@testset "element ops" begin
    frames = framed_signal(x,fs,0.09,0.01)
    frame = extract_frame(frames,11)
    mag = magspec(frame,fs)
    msp = specgram(frames)

    #Each op maps the components and leaves everything else untouched
    lmag = log(mag)
    l10mag = log10(mag)
    dbmag = amp2db(mag)
    pdbmag = pow2db(mag)
    @test lmag.components == log.(mag.components)
    @test l10mag.components == log10.(mag.components)
    @test dbmag.components == amp2db.(mag.components)
    @test pdbmag.components == pow2db.(mag.components)
    @test lmag.frqs == mag.frqs
    @test lmag.signal == mag.signal
    @test lmag.title == mag.title

    lmsp = log(msp)
    l10msp = log10(msp)
    dbmsp = amp2db(msp)
    pdbmsp = pow2db(msp)
    @test lmsp.components == log.(msp.components)
    @test l10msp.components == log10.(msp.components)
    @test dbmsp.components == amp2db.(msp.components)
    @test pdbmsp.components == pow2db.(msp.components)
    @test lmsp.frqs == msp.frqs
    @test lmsp.time == msp.time
    @test lmsp.frames === msp.frames
end
