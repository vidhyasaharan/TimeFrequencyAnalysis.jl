#Spectral analyses on signals with known content: pure tones must produce their peaks in
#the correct bins with the predicted magnitudes, the comp() forms must return the same
#numbers as the container forms, and the frequency axes must map bins to Hz correctly.

@testset "magspec" begin
    #100 samples of a full-scale 25 Hz cosine at fs = 100 Hz with a rectangular window:
    #the one-sided spectrum has 100/2 + 1 = 51 bins at 1 Hz per bin (DC in bin 1), the
    #peak must land in bin 26 (= 25 Hz), and an on-bin unit cosine has magnitude N/2 = 50
    t = 0:0.01:0.99
    frq = 25
    mfs = 100.0
    xx = cos.(2*pi*frq*t)
    spec = magspec(xx,mfs;wtype = "rect")
    mspec1, frqs1 = magspec(comp(), xx, mfs; wtype = "rect")
    mspec = spec.components
    @test mspec1 == mspec #comp() form returns the same numbers as the container form
    @test frqs1 == spec.frqs
    mmag,mfrq = findmax(mspec)
    @test length(mspec) == 51
    @test mfrq == frq + 1 #peak at the tone frequency (1 Hz per bin, DC in bin 1)
    @test round(mmag) == 50.0

    #An ndft below the signal length is raised to the signal length, grid included
    mspec2, frqs2 = magspec(comp(), xx, mfs; ndft = 10, wtype = "rect")
    @test length(mspec2) == length(frqs2)
    @test length(mspec2) == 51
end

@testset "cexp" begin
    #cexp(f, fs, N) must be the unit-norm complex exponential probe e^(-i2πfn/fs)/√N:
    #after undoing the 1/√N scaling, the real part is cos and the imaginary part -sin at
    #every sample, for on-bin and off-bin frequencies alike
    cfs = 8000
    N = 16000
    cs(fr,i) = cos(2π*(fr/cfs)*i)
    sn(fr,i) = sin(2π*(fr/cfs)*i)

    for f in [100, 700, 3401]
        ce = sqrt(N)*TimeFrequencyAnalysis.cexp(f,cfs,N)
        @test ce⋅ce ≈ N #unit norm before scaling
        for i=1:N
            @test isapprox(real(ce[i]),cs(f,i);atol = 1e-9)
            @test isapprox(imag(ce[i]),-sn(f,i);atol = 1e-9)
        end
    end

    @test eltype(TimeFrequencyAnalysis.cexp(Float32, 100, cfs, 16)) == Complex{Float32}
end

@testset "cexp_proj_matrix" begin
    #The projection matrix must stack one cexp probe vector per requested frequency as
    #its rows (frequencies × samples), so periodogram is a single matrix-vector product
    cfs = 16000
    frqs = [10, 50, 100, 500, 1000, 5000]
    N = 1600
    proj = TimeFrequencyAnalysis.cexp_proj_matrix(frqs,cfs,N)
    @test typeof(proj) <: Matrix{Complex{Float64}}
    @test size(proj,1) == length(frqs)
    @test size(proj,2) == N
    for i = 1:length(frqs)
        ce = TimeFrequencyAnalysis.cexp(frqs[i],cfs,N)
        @test proj[i,:] == ce
    end
end

@testset "periodogram" begin
    #A 1024 Hz tone analysed on an explicit frequency grid that includes 1024 Hz: the
    #periodogram must peak exactly at that grid entry, with one value per grid frequency
    pfs = 8000
    n = 1:160
    frq = 1024
    xx = cos.(2*pi*(frq/pfs)*n)

    tfrqs = [100, 500, 705, 1024, 1800, 2100, 3401, 3700]
    spec = periodogram(xx,pfs,tfrqs)
    @test spec isa spectrum{Float64,Float64}
    pgram  = spec.components
    @test pgram == periodogram(comp(), xx,pfs,tfrqs)
    @test length(pgram) == length(tfrqs)
    mmag, mfindx = findmax(pgram)
    @test tfrqs[mfindx] == frq #the periodogram peaks at the tone frequency

    #The fmin/fmax form evaluates on the default logarithmic grid over that range and
    #must still peak at the tone (1024 = 2^10 lies exactly on the log2 grid from 8 Hz)
    spec = periodogram(xx,pfs;fmin=8,fmax=4000)
    @test spec isa spectrum{Float64,Float64}
    pgram = spec.components
    @test pgram == periodogram(comp(),xx,pfs,fmin=8,fmax=4000)
    mmag, mfindx = findmax(pgram)
    frqs = logfreq_array(;fmin=8,fmax=4000)
    @test length(pgram) == length(frqs)
    @test frqs[mfindx] == frq
end

@testset "specgram" begin
    #On the shared broadband fixture: the spectrogram must be a real positive matrix
    #(floored at eps to keep log/dB well defined) with one frequency per row on an axis
    #running from DC to (at most) Nyquist
    tf = specgram(x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    msp = tf.components
    @test msp == specgram(comp(), x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    @test typeof(msp) == Matrix{Float64}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test length(tf.frqs) == size(msp,1) #one frequency index per spectrogram row
    @test tf.frqs[1] == 0.0
    @test (tf.frqs[end] < fs/2) || (tf.frqs[end] ≈ fs/2) #Nyquist bin up to rounding in fs/nfft
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()

    #Two one-second tones back to back with one-second non-overlapping frames: each tone
    #must dominate its own column, peaking in the bin of its frequency (1 Hz per bin,
    #DC in row 1), and the frequency axis must map those rows back to Hz
    t = 0:0.01:0.99
    frq1 = 25
    frq2 = 30
    xx = [cos.(2*pi*frq1*t);cos.(2*pi*frq2*t)]
    tf = specgram(xx,100.0,frame_dur = 1.0, frame_shift_dur = 1.0)
    msp = tf.components
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()
    mmag1,mfrq1 = findmax(msp[:,1])
    mmag2,mfrq2 = findmax(msp[:,2])
    @test mfrq1 == frq1 + 1
    @test mfrq2 == frq2 + 1
    @test length(tf.frqs) == size(msp,1)
    @test tf.frqs[mfrq1] ≈ frq1 #frequency axis maps spectrogram rows to Hz
    @test tf.frqs[mfrq2] ≈ frq2
end


#welch_psd is a power estimator, not a magnitude one: the two scaling conventions are exact
#statements about what a bin means, and both are asserted directly rather than by comparison
#with another function.
@testset "welch_psd scaling conventions" begin
    wfs = 1000.0
    N = 200_000
    t = (0:N-1)./wfs

    #:spectrum - a sinusoid of amplitude A reads A^2/2 at its bin, INDEPENDENT of frame length.
    #This is the convention for finding discrete lines: the number a line reports does not move
    #when the analysis resolution changes.
    A = 1.7
    tone = A.*cos.(2π*50 .*t)
    for fd in (0.5, 1.0, 2.0)
        sp = welch_psd(tone, wfs; frame_dur = fd, scaling = :spectrum)
        @test maximum(sp.components) ≈ A^2/2 rtol = 1e-4
        @test sp.frqs[argmax(sp.components)] ≈ 50
    end

    #:density - Parseval: summing the density over the frequency grid recovers the variance,
    #again independently of frame length. This is the convention for noise floors.
    σ = 2.0
    noise = σ.*randn(N)
    for fd in (0.5, 1.0, 2.0)
        sp = welch_psd(noise, wfs; frame_dur = fd, scaling = :density)
        Δf = sp.frqs[2] - sp.frqs[1]
        @test sum(sp.components)*Δf ≈ σ^2 rtol = 0.02
    end

    #A white noise density is flat at 2σ²/fs across the band
    sp = welch_psd(noise, wfs; frame_dur = 1.0, scaling = :density)
    interior = sp.components[10:end-10]
    @test sum(interior)/length(interior) ≈ 2σ^2/wfs rtol = 0.05

    #Averaging is what makes the estimate reliable: the scatter of a noise density falls as
    #1/sqrt(F) while its mean stays put. A single frame has ~100% standard error however long
    #it is, which is the whole reason this function exists rather than magspec.
    scatter = Float64[]
    for nfr in (1, 16, 256)
        dur = nfr == 1 ? N/wfs : (N/wfs)/((nfr+1)/2)   #50% overlap => F ≈ 2*total/frame - 1
        sp = welch_psd(noise, wfs; frame_dur = dur, scaling = :density)
        b = sp.components[10:end-10]
        m = sum(b)/length(b)
        push!(scatter, sqrt(sum(abs2, b .- m)/length(b))/m)
        @test m ≈ 2σ^2/wfs rtol = 0.1                  #the mean does not move
    end
    @test scatter[1] > 0.5                             #one frame: standard error of order 1
    @test scatter[3] < scatter[2] < scatter[1]         #and it falls monotonically with F
    @test scatter[3] < 0.2
end


#The four dc modes are two independent operations. Dividing is exact and window independent;
#subtracting is not a DC-bin-only operation once a tapered window is involved, and suppressing
#that leakage is the point of it.
@testset "welch_psd dc handling" begin
    wfs = 100.0
    N = 6000
    t = (0:N-1)./wfs
    μ = 2.5
    env = μ.*(1 .+ 0.4.*cos.(2π*3 .*t))                #strictly positive, constant frame mean
    kw = (frame_dur = 4.0, scaling = :spectrum)

    Pk = welch_psd(comp(), env, wfs; dc = :keep,   kw...)[1]
    Pr = welch_psd(comp(), env, wfs; dc = :remove, kw...)[1]
    Pd = welch_psd(comp(), env, wfs; dc = :divide, kw...)[1]
    Pi = welch_psd(comp(), env, wfs; dc = :index,  kw...)[1]

    #Division rescales every bin by exactly 1/mu^2, whatever the window
    for wt in ("hanning", "rect", "hamming")
        a = welch_psd(comp(), env, wfs; dc = :keep,   wtype = wt, kw...)[1]
        b = welch_psd(comp(), env, wfs; dc = :remove, wtype = wt, kw...)[1]
        c = welch_psd(comp(), env, wfs; dc = :divide, wtype = wt, kw...)[1]
        d = welch_psd(comp(), env, wfs; dc = :index,  wtype = wt, kw...)[1]
        @test c ≈ a./μ^2
        @test d ≈ b./μ^2
    end

    #Subtracting the mean reaches only the DC bin for a rectangular window...
    ak = welch_psd(comp(), env, wfs; dc = :keep,   wtype = "rect", kw...)[1]
    ar = welch_psd(comp(), env, wfs; dc = :remove, wtype = "rect", kw...)[1]
    @test ak[2:end] ≈ ar[2:end]
    @test ak[1] > 1e6*ar[1]

    #...but a tapered window spreads it over the mainlobe, which is exactly why the mean is
    #removed before windowing rather than left for the window to contain. The bin next to DC
    #carries the mean's leakage under :keep and essentially nothing under :remove.
    @test Pk[2] > 1e6*Pr[2]
    @test Pr[2] < 1e-6

    #The modulation line itself is untouched by the dc rule: same bin, and the dimensionless
    #modes report it a factor mu^2 lower
    kline = argmax(Pk[3:end]) + 2
    @test argmax(Pr[3:end]) + 2 == kline
    @test argmax(Pi[3:end]) + 2 == kline
    @test Pk[kline] ≈ Pi[kline]*μ^2 rtol = 1e-3
    #a 40% modulation of a mean-mu envelope reads (0.4)^2/2 as a modulation index
    @test Pi[kline] ≈ 0.4^2/2 rtol = 1e-3

    #The dimensionless modes are invariant to overall level; the physical ones scale with it
    Pi10 = welch_psd(comp(), 10 .*env, wfs; dc = :index, kw...)[1]
    Pk10 = welch_psd(comp(), 10 .*env, wfs; dc = :keep,  kw...)[1]
    @test Pi10[kline] ≈ Pi[kline] rtol = 1e-6
    @test Pk10[kline] ≈ 100*Pk[kline] rtol = 1e-6
end


@testset "welch_psd containers and argument checking" begin
    wfs = 1000.0
    xw = randn(4000)
    sw = signal(xw, wfs)

    #Container form carries the signal, the axis and a title naming the scaling
    sp = welch_psd(sw; frame_dur = 0.5)
    @test sp isa spectrum{Float64,Float64}
    @test length(sp.components) == length(sp.frqs)
    @test sp.frqs[1] == 0
    @test sp.title == "Welch Power Spectrum"
    @test welch_psd(sw; frame_dur = 0.5, scaling = :density).title == "Welch Power Spectral Density"

    #comp() forms: bare vector from a framed_signal (as specgram), (P, frqs) from the signal
    #and array forms (as magspec)
    frames = framed_signal(sw, 0.5, 0.25)
    P = welch_psd(comp(), frames)
    @test P isa Vector{Float64}
    @test P == sp.components
    Pt, ft = welch_psd(comp(), sw; frame_dur = 0.5)
    @test Pt == sp.components && ft == sp.frqs
    @test welch_psd(comp(), xw, wfs; frame_dur = 0.5)[1] == sp.components
    @test welch_psd(xw, wfs; frame_dur = 0.5).components == sp.components
    @test welch_psd(frames).components == sp.components

    #Only whole frames are averaged, so a zero padded tail frame cannot drag the estimate down
    @test frames.num_signal_frames < frames.num_frames

    #ndft is raised to the frame length, never below it
    @test length(welch_psd(comp(), framed_signal(sw, 0.5, 0.25); ndft = 8)) == div(500,2) + 1

    #Precision follows the signal
    @test welch_psd(signal(randn(Float32, 4000), 1000f0); frame_dur = 0.5) isa spectrum{Float32,Float32}

    @test_throws ErrorException welch_psd(sw; frame_dur = 0.5, dc = :normalise)
    @test_throws ErrorException welch_psd(sw; frame_dur = 0.5, scaling = :power)
    @test_throws ErrorException welch_psd(signal(randn(100), wfs); frame_dur = 0.5)  #shorter than a frame
    #:divide and :index divide by the frame mean, so they refuse a frame whose mean is
    #negligible against its own RMS. A tone spanning whole periods in every frame is exactly
    #that case: 25 periods per 0.5 s frame at 50 Hz leaves a mean at the rounding level.
    centred = signal(cos.(2π*50 .*(0:3999)./wfs), wfs)
    @test_throws ErrorException welch_psd(centred; frame_dur = 0.5, dc = :index)
    @test_throws ErrorException welch_psd(centred; frame_dur = 0.5, dc = :divide)
    #An all-zero frame is not an error: it simply contributes no power under any rule
    @test all(welch_psd(comp(), signal(zeros(4000), wfs); frame_dur = 0.5, dc = :index)[1] .== 0)
end
