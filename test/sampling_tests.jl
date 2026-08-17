#Duration/sample-count conversion and resampling: time2nsamples must round to the
#nearest whole sample, and resample must change length and rate consistently while
#preserving the signal content.

@testset "time2nsamples" begin
    #Exact products convert exactly (0.02 s at 16 kHz is 320 samples; 1 s at 8 kHz is
    #8000); non-integer products round to the nearest sample; the result is always an Int
    @test time2nsamples(0.02, 16000) == 320
    @test time2nsamples(1, 8000.0) == 8000
    @test time2nsamples(0.0075, 1000) == 8 #7.5 samples rounds to nearest
    @test typeof(time2nsamples(0.01, 44100.0)) == Int
end


@testset "resample" begin
    #Doubling the sampling rate of a 1 kHz tone must double the number of samples (up to
    #±1 from the polyphase filter edges) and report the new rate
    rfs = 8000
    n = 1:2*rfs
    frq = 1000
    xx = cos.(2*pi*(frq/rfs)*n)

    s = signal(xx,rfs)
    rs = resample(s,2*rfs)
    @test rs isa signal{Float64}
    @test abs(length(rs.x) - length(s.x)*2) <= 1
    @test rs.fs == s.fs*2

    #The dominant frequency must be preserved by resampling: the magnitude-spectrum peak
    #of the resampled signal must be the grid frequency closest to the original peak
    #(closer than either neighbouring bin)
    mag = magspec(s)
    rmag = magspec(rs)

    ind = argmax(mag.components)
    rind = argmax(rmag.components)

    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind-1])
    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind+1])
end


#The rate decides the algorithm: an exact small fraction takes the rational polyphase path (one
#branch per output sample), anything else falls back to the arbitrary-rate resampler.
@testset "resample rate resolution" begin
    #Integer decimation and interpolation are exact, whichever way the rates are written
    @test TimeFrequencyAnalysis._resample_rate(3125, 312.5) === 1//10
    @test TimeFrequencyAnalysis._resample_rate(8000, 16000) === 2//1
    @test TimeFrequencyAnalysis._resample_rate(44100, 48000) === 160//147

    #An incommensurable pair cannot be written as a small fraction and must not be forced into
    #one - it falls back to Float64 rather than silently asking for thousands of phases
    @test TimeFrequencyAnalysis._resample_rate(1000, 1000π) isa Float64

    @test_throws ErrorException TimeFrequencyAnalysis._resample_rate(0, 100)
    @test_throws ErrorException TimeFrequencyAnalysis._resample_rate(100, -1)
end


#The array methods exist so a whole subband matrix can be decimated in one call (see
#power_envelope), and the precision fix applies to every form: DSP designs its taps in Float64
#and its output promotes to match them, so without converting the taps a Float32 input would
#come back Float64.
@testset "resample arrays and precision" begin
    afs, afs_new = 3125.0, 312.5
    N = 2000

    #Each row of the matrix form is the vector form of that row, and the vector form of a
    #signal's samples is the signal form's samples
    X = randn(3, N)
    R = resample(X, afs, afs_new; dims = 2)
    @test size(R) == (3, ceil(Int, N*afs_new/afs))
    for k in axes(X,1)
        @test R[k,:] == resample(X[k,:], afs, afs_new)
    end
    @test resample(X[1,:], afs, afs_new) == resample(signal(X[1,:], afs), afs_new).x

    #dims selects the axis, and resampling columns is resampling the transpose's rows
    Xc = permutedims(X)
    @test resample(Xc, afs, afs_new; dims = 1) ≈ permutedims(R)
    @test_throws ErrorException resample(X, afs, afs_new; dims = 3)

    #The anti-alias filter has unit DC gain, so decimation preserves the level exactly. This is
    #the invariant power_envelope relies on: a decimated envelope has the same mean as the
    #envelope it came from. Away from the edges the constant is reproduced to machine precision;
    #the transient at each end is the filter's half-length expressed at the decimated rate, so
    #derive the margin from the taps rather than guessing at it.
    ntaps = length(TimeFrequencyAnalysis._resample_taps(Float64, 1//10, 1, 60))
    margin = ceil(Int, (ntaps - 1)/2*(afs_new/afs)) + 1
    c = resample(fill(2.5, 5000), afs, afs_new)
    @test all(isapprox.(c[margin+1:end-margin], 2.5; atol = 1e-12))
    @test sum(c[margin+1:end-margin])/length(c[margin+1:end-margin]) ≈ 2.5
    #and the transient really is confined to that margin - one sample inside it, it is not
    @test !all(isapprox.(c[margin-9:end-margin], 2.5; atol = 1e-12))

    #Precision propagates through every form
    X32 = randn(Float32, 3, N)
    @test eltype(resample(X32, afs, afs_new; dims = 2)) == Float32
    @test eltype(resample(X32[1,:], afs, afs_new)) == Float32
    @test resample(signal(X32[1,:], Float32(afs)), afs_new) isa signal{Float32}
    @test resample(signal(randn(N), afs), afs_new) isa signal{Float64}

    #Float32 and Float64 agree to single precision
    @test resample(X32, afs, afs_new; dims = 2) ≈ Float32.(resample(Float64.(X32), afs, afs_new; dims = 2)) rtol = 1e-4

    #A tone survives decimation: 40 Hz is well inside the 156.25 Hz Nyquist of the new rate
    t = (0:N-1)./afs
    tone = cos.(2π*40 .*t)
    d = resample(tone, afs, afs_new)
    m, f = magspec(comp(), d[20:end-20], afs_new)
    @test abs(f[argmax(m)] - 40) < 1.0
end
