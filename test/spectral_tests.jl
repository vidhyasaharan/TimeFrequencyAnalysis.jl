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
