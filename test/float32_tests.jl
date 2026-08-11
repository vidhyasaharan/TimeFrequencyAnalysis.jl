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

@testset "correlations" begin
    #Correlation outputs must stay in the input precision, for the sequence functions
    #and for both the single-point and batch (matrix) forms of nacf
    @test eltype(xcorr(randn(Float32,100), randn(Float32,10))) == Float32
    @test eltype(acorr(randn(Float32,1000), 10)) == Float32
    s32 = signal(randn(Float32, 2000), 8000)
    @test nacf(s32, 1, 40; win_size = 200) isa Float32
    @test eltype(nacf(s32; min_lag = 1, max_lag = 20, win_size = 200, win_shift = 100)) == Float32
end
