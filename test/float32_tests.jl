#Verify that the whole chain is generic over the floating point precision of the input:
#a Float32 signal must produce Float32 containers and results at every step, with no
#silent promotion to Float64 anywhere in between.
@testset "objects" begin
    x32 = randn(Float32, 16000)
    s32 = signal(x32, 16000)
    @test s32 isa signal{Float32}
    @test typeof(s32.fs) == Float32

    frames32 = framed_signal(s32)
    @test frames32 isa framed_signal{Float32}
    @test eltype(extract_frame(frames32,1)) == Float32
    @test eltype(frame_energy(frames32)) == Float32
end

@testset "framing and windows" begin
    x32 = randn(Float32, 1600)
    @test eltype(window(Float32, 32)) == Float32
    @test window(Float32, 32) ≈ window(32)
    @test eltype(enframe(x32, 160, 80)) == Float32
end

@testset "spectral analyses" begin
    x32 = randn(Float32, 16000)
    s32 = signal(x32, 16000)
    frames32 = framed_signal(s32)

    @test eltype(dft(randn(Float32,256))) == Complex{Float32}

    ms = magspec(s32)
    @test ms isa spectrum{Float32,Float32}
    @test eltype(ms.frqs) == Float32

    tf = specgram(frames32)
    @test tf isa timefreq{Float32,Float32}
    @test length(tf.frqs) == size(tf.components,1)
    @test eltype(tf.time) == Float32
    @test amp2db(tf) isa timefreq{Float32,Float32}

    pg = periodogram(comp(), frames32)
    @test eltype(pg) == Float32

    #Float32 and Float64 paths must agree numerically
    x64 = convert(Vector{Float64}, x32)
    tf64 = specgram(framed_signal(x64, 16000.0))
    @test isapprox(tf.components, tf64.components; rtol = 1e-3)
end

@testset "gammatone" begin
    #The gammatone chain computes in the signal's precision even through a Float64 bank,
    #and the Float32 and Float64 paths agree numerically
    x32 = randn(Float32, 8000)
    s32 = signal(x32, 16000)
    fb = gammatone_filterbank(16000.0)
    tf = gammatone_analysis(s32, fb)
    @test tf isa timefreq{Float32, ComplexF32}
    @test eltype(tf.time) == Float32
    cg = gammatone_cochleagram(s32, fb)
    @test cg isa timefreq{Float32, Float32}
    @test amp2db(cg) isa timefreq{Float32, Float32}

    cg64 = gammatone_cochleagram(signal(convert(Vector{Float64}, x32), 16000.0), fb)
    @test isapprox(cg.components, cg64.components; rtol = 1e-3)
end

@testset "resampling" begin
    #DSP designs its resampling taps in Float64 and its output promotes to match them, so
    #without converting the taps to the data's precision every resampled Float32 signal would
    #come back Float64. This is the regression test for that conversion, on all three forms.
    x32 = randn(Float32, 4000)
    rs = resample(signal(x32, 8000), 4000)
    @test rs isa signal{Float32}
    @test typeof(rs.fs) == Float32
    @test eltype(resample(x32, 8000, 4000)) == Float32

    X32 = randn(Float32, 3, 4000)
    @test eltype(resample(X32, 8000, 4000; dims = 2)) == Float32
    @test eltype(resample(permutedims(X32), 8000, 4000; dims = 1)) == Float32

    #Float32 and Float64 paths must agree numerically
    @test isapprox(rs.x, convert(Vector{Float64}, x32) |> x -> resample(signal(x, 8000.0), 4000).x; rtol = 1e-3)
end

@testset "modulation" begin
    #The whole modulation chain in Float32: the filterbank is designed at Float32, the envelope
    #and its decimation stay there, and the modfreq that comes out carries Float32 axes. This is
    #the end-to-end guard - the intermediate stages each have their own precision tests, but only
    #this one runs the full signal -> subbands -> envelope -> decimate -> Welch path.
    mfs32 = 3125f0
    t = (0:round(Int, 20*3125)-1)./3125
    x32 = Float32.((1 .+ 0.5cos.(2π*7 .*t)).*cos.(2π*200 .*t))
    s32 = signal(x32, mfs32)
    fb32 = gammatone_filterbank(mfs32; fmin = 100, fmax = 400, base_frq = 200,
                                filters_per_erb = 2, bw_scale = 2)
    @test fb32 isa gammatone_filterbank{Float32}
    @test default_align(fb32) isa Float32

    env = power_envelope(s32, fb32)
    @test env isa timefreq{Float32,Float32}
    @test eltype(env.time) == Float32
    @test eltype(power_envelope(comp(), s32, fb32)) == Float32

    ms = modulation_spectrum(s32, fb32; frame_dur = 4.0)
    @test ms isa modfreq{Float32,Float32}
    @test eltype(ms.fcs) == Float32
    @test eltype(ms.mod_frqs) == Float32
    @test eltype(ms.bws) == Float32
    @test typeof(ms.fs_env) == Float32
    @test typeof(ms.frame_dur) == Float32
    @test pow2db(ms) isa modfreq{Float32,Float32}

    #The Float32 and Float64 paths must agree numerically, and must find the same cell
    fb64 = gammatone_filterbank(3125.0; fmin = 100, fmax = 400, base_frq = 200,
                                filters_per_erb = 2, bw_scale = 2)
    ms64 = modulation_spectrum(signal(convert(Vector{Float64}, x32), 3125.0), fb64; frame_dur = 4.0)
    @test isapprox(ms.components, ms64.components; rtol = 1e-2)
    @test Tuple(argmax(ms.components)) == Tuple(argmax(ms64.components))

    #welch_psd on its own, for each dc mode, since :index divides by a Float32 mean
    env32 = power_envelope(comp(), s32, fb32)[1,:]
    for d in (:keep, :remove, :divide, :index)
        @test welch_psd(env32, ms.fs_env; frame_dur = 4.0, dc = d) isa spectrum{Float32,Float32}
    end
end

@testset "correlations" begin
    #Correlation outputs must stay in the input precision, for the sequence functions
    #and for both the single-point and batch (matrix) forms of nacf
    @test eltype(xcorr(randn(Float32,100), randn(Float32,10))) == Float32
    @test eltype(acorr(randn(Float32,1000), 10)) == Float32
    s32 = signal(randn(Float32, 2000), 8000)
    @test nacf(s32, 1, 40; win_size = 200) isa Float32
    @test eltype(nacf(s32; min_lag = 1, max_lag = 20, win_size = 200, win_shift = 100)) == Float32
end
