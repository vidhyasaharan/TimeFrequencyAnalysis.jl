#Synthetic signal generators (the package's test fixtures): correct output length for
#the sample-count and duration call forms, requested element type, and the documented
#structure of each signal (fresh noise per call, impulse positions).

@testset "white_noise" begin
    #Length and element type follow the arguments; the duration form spans dur*fs samples
    wn = white_noise(1000)
    @test typeof(wn) == Vector{Float64}
    @test length(wn) == 1000
    @test typeof(white_noise(Float32, 10)) == Vector{Float32}
    @test length(white_noise(0.5, 1000)) == 500 #duration form: 0.5 s at 1 kHz
    #Each call draws from a freshly entropy-seeded generator, so two calls must differ
    @test white_noise(1000) != wn
end

@testset "ar_process" begin
    #White noise filtered through a stable all-pole model: output length follows the
    #sample-count and duration call forms
    a = [1.0, 0.0, 0.0, 0.8] #stable AR(3)
    ap = ar_process(a, 5000)
    @test typeof(ap) == Vector{Float64}
    @test length(ap) == 5000
    @test length(ar_process(a, 0.5, 1000)) == 500
end

@testset "impulse_train" begin
    #Period 10 over 95 samples: unit impulses at samples 1, 11, ..., 91 (10 in total),
    #zeros elsewhere
    it = impulse_train(10, 95)
    @test typeof(it) == Vector{Float64}
    @test length(it) == 95
    @test sum(it) == 10 #impulses at samples 1, 11, ..., 91
    @test it[1] == 1
    @test it[2] == 0
    @test it[11] == 1
    @test typeof(impulse_train(Float32, 10, 95)) == Vector{Float32}

    #Frequency/duration form: 100 Hz for 1 s at 8 kHz is an 80-sample period, giving
    #8000 samples containing exactly 100 impulses
    it2 = impulse_train(100, 1.0, 8000)
    @test length(it2) == 8000
    @test sum(it2) == 100
end
