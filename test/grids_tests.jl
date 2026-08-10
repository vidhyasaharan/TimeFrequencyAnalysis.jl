@testset "logfreq_array" begin
    frqs = logfreq_array(;fmin = 2, fmax = 8192, frq_per_octave = 100)
    lfrqs = log2.(frqs)
    @test lfrqs[1] == 1.0
    @test lfrqs[end] == 13.0
    @test length(frqs) == 1201
    for i=1:length(frqs)
        @test frqs[i] ≈ exp2(0.99 + 0.01*i)
    end
end

@testset "linfreq_array" begin
    frqs = linfreq_array(;fmin = 0, fmax = 100, nfrqs = 11)
    @test frqs[1] == 0.0
    @test frqs[end] == 100.0
    for i=2:length(frqs)
        @test frqs[i] == (i-1)*10
    end
end

@testset "findclosest" begin
    data = 0.0:0.1:10.0
    rin = 7
    xc = (rin-1)*0.1
    cin = TimeFrequencyAnalysis.findclosest(xc,data)
    @test cin == rin

    #Mixed input types are supported
    @test TimeFrequencyAnalysis.findclosest(3, collect(0.0:1.0:10.0)) == 4
    @test TimeFrequencyAnalysis.findclosest(2.4f0, collect(0.0:1.0:10.0)) == 3
end

@testset "frqindex" begin
    frqs = 0.0:10.0:1000.0
    f = 30.0
    fin = frqindex(f,frqs)
    @test typeof(fin) == Int
    @test fin == 4

    fv = [21.0, 41.5, 59.0]
    fvin = frqindex(fv,frqs)
    @test typeof(fvin) == Vector{Int}
    @test length(fvin) == length(fv)
    @test fvin[1] == 3
    @test fvin[2] == 5
    @test fvin[3] == 7
end
