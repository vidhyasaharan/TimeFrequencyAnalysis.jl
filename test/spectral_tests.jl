@testset "magspec" begin
    t = 0:0.01:0.99
    frq = 25
    mfs = 100.0
    xx = cos.(2*pi*frq*t)
    spec = magspec(xx,mfs;wtype = "rect")
    mspec1, frqs1 = magspec(comp(), xx, mfs; wtype = "rect")
    mspec = spec.components
    @test mspec1 == mspec
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
    cfs = 16000
    frqs = [10, 50, 100, 500, 1000, 5000]
    N = 1600
    proj = TimeFrequencyAnalysis.cexp_proj_matrix(frqs,cfs,N)
    @test typeof(proj) <: Matrix{Complex{Float}}
    @test size(proj,1) == length(frqs)
    @test size(proj,2) == N
    for i = 1:length(frqs)
        ce = TimeFrequencyAnalysis.cexp(frqs[i],cfs,N)
        @test proj[i,:] == ce
    end
end

@testset "periodogram" begin
    pfs = 8000
    n = 1:160
    frq = 1024
    xx = cos.(2*pi*(frq/pfs)*n)

    tfrqs = [100, 500, 705, 1024, 1800, 2100, 3401, 3700]
    spec = periodogram(xx,pfs,tfrqs)
    @test spec isa spectrum{Float,Float}
    pgram  = spec.components
    @test pgram == periodogram(comp(), xx,pfs,tfrqs)
    @test length(pgram) == length(tfrqs)
    mmag, mfindx = findmax(pgram)
    @test tfrqs[mfindx] == frq #the periodogram peaks at the tone frequency

    spec = periodogram(xx,pfs;fmin=8,fmax=4000)
    @test spec isa spectrum{Float,Float}
    pgram = spec.components
    @test pgram == periodogram(comp(),xx,pfs,fmin=8,fmax=4000)
    mmag, mfindx = findmax(pgram)
    frqs = logfreq_array(;fmin=8,fmax=4000)
    @test length(pgram) == length(frqs)
    @test frqs[mfindx] == frq
end

@testset "specgram" begin
    tf = specgram(x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    msp = tf.components
    @test msp == specgram(comp(), x, fs;frame_dur = 0.03,frame_shift_dur=0.01, wtype = "hamming")
    @test typeof(msp) == Matrix{Float}
    @test size(msp,1) > 0
    @test size(msp,2) > 0
    @test length(tf.frqs) == size(msp,1) #one frequency index per spectrogram row
    @test tf.frqs[1] == 0.0
    @test (tf.frqs[end] < fs/2) || (tf.frqs[end] ≈ fs/2) #Nyquist bin up to rounding in fs/nfft
    @test maximum(isa.(msp,Complex))==false
    @test minimum(msp) >= eps()

    #Two one-second tones back to back with one-second frames land in separate columns
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
