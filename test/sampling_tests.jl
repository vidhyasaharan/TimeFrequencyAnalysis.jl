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
