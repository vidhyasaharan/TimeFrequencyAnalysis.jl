@testset "white_noise" begin
    wn = white_noise(1000)
    @test typeof(wn) == Vector{Float}
    @test length(wn) == 1000
    @test typeof(white_noise(Float32, 10)) == Vector{Float32}
    @test length(white_noise(0.5, 1000)) == 500 #duration form
    @test white_noise(1000) != wn #each call draws fresh noise
end

@testset "ar_process" begin
    a = [1.0, 0.0, 0.0, 0.8] #stable AR(3)
    ap = ar_process(a, 5000)
    @test typeof(ap) == Vector{Float}
    @test length(ap) == 5000
    @test length(ar_process(a, 0.5, 1000)) == 500
end

@testset "impulse_train" begin
    it = impulse_train(10, 95)
    @test typeof(it) == Vector{Float}
    @test length(it) == 95
    @test sum(it) == 10 #impulses at samples 1, 11, ..., 91
    @test it[1] == 1
    @test it[2] == 0
    @test it[11] == 1
    @test typeof(impulse_train(Float32, 10, 95)) == Vector{Float32}

    it2 = impulse_train(100, 1.0, 8000) #f₀ = 100 Hz for 1 s at 8 kHz
    @test length(it2) == 8000
    @test sum(it2) == 100
end
