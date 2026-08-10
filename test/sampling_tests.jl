@testset "time2nsamples" begin
    @test time2nsamples(0.02, 16000) == 320
    @test time2nsamples(1, 8000.0) == 8000
    @test time2nsamples(0.0075, 1000) == 8 #rounds to nearest
    @test typeof(time2nsamples(0.01, 44100.0)) == Int
end


@testset "resample" begin
    rfs = 8000
    n = 1:2*rfs
    frq = 1000
    xx = cos.(2*pi*(frq/rfs)*n)

    s = waveform(xx,rfs)
    rs = resample(s,2*rfs)
    @test rs isa waveform{Float}
    @test abs(length(rs.x) - length(s.x)*2) <= 1
    @test rs.fs == s.fs*2

    #The dominant frequency must be preserved by resampling
    mag = magspec(s)
    rmag = magspec(rs)

    ind = argmax(mag.components)
    rind = argmax(rmag.components)

    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind-1])
    @test abs(mag.frqs[ind] - rmag.frqs[rind]) < abs(mag.frqs[ind] - rmag.frqs[rind+1])
end
