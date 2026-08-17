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


#The same element-wise contract for modfreq: components transformed, every axis and every piece
#of provenance carried through untouched.
@testset "modfreq element ops" begin
    C = [1.0 10.0; 100.0 1000.0]
    mf = modfreq(C, [80.0, 160.0], [0.0, 5.0], [66.0, 83.0], 312.5, 10.0, 11, "Modulation Spectrum")

    for (f, expected) in ((log, map(log, C)), (log10, map(log10, C)),
                          (amp2db, map(amp2db, C)), (pow2db, map(pow2db, C)))
        out = f(mf)
        @test out isa modfreq{Float64,Float64}
        @test out.components == expected
        @test out.fcs == mf.fcs
        @test out.mod_frqs == mf.mod_frqs
        @test out.bws == mf.bws
        @test out.fs_env == mf.fs_env
        @test out.frame_dur == mf.frame_dur
        @test out.nframes == mf.nframes
        @test out.title == mf.title
    end

    #pow2db is 10log10 and amp2db 20log10, so one is exactly half the other on the same data
    @test pow2db(mf).components ≈ amp2db(mf).components ./ 2
    @test eltype(pow2db(modfreq(Float32.(C), [80.0, 160.0], [0.0, 5.0], [66.0, 83.0], 312.5, 10.0, 11)).components) == Float32
end
