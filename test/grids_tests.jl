#Frequency grids and nearest-value indexing: grids must span exactly [fmin, fmax] with
#the documented spacing (per-octave for logarithmic, uniform for linear), and the index
#lookups must return the position of the closest grid value.

@testset "logfreq_array" begin
    #2 Hz to 8192 Hz is octaves 1 to 13; at 100 points per octave that is 12*100 + 1 =
    #1201 points, and point i must sit at 2^(0.99 + 0.01i) (uniform spacing in log2)
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
    #0 to 100 Hz in 11 points: end points included, uniform 10 Hz spacing
    frqs = linfreq_array(;fmin = 0, fmax = 100, nfrqs = 11)
    @test frqs[1] == 0.0
    @test frqs[end] == 100.0
    for i=2:length(frqs)
        @test frqs[i] == (i-1)*10
    end
end

@testset "findclosest" begin
    #findclosest must return the index of the grid value nearest to the query; a query
    #that lies exactly on the grid returns its own position
    data = 0.0:0.1:10.0
    rin = 7
    xc = (rin-1)*0.1
    cin = TimeFrequencyAnalysis.findclosest(xc,data)
    @test cin == rin

    #Mixed input types are supported: 3 (Int) is exactly grid value 3.0 at index 4;
    #2.4f0 (Float32) is nearest to 2.0 at index 3
    @test TimeFrequencyAnalysis.findclosest(3, collect(0.0:1.0:10.0)) == 4
    @test TimeFrequencyAnalysis.findclosest(2.4f0, collect(0.0:1.0:10.0)) == 3
end

@testset "frqindex" begin
    #frqindex maps frequencies to nearest-grid-value indices on the grid 0:10:1000
    #(index 1 is 0 Hz): 30 Hz sits exactly at index 4
    frqs = 0.0:10.0:1000.0
    f = 30.0
    fin = frqindex(f,frqs)
    @test typeof(fin) == Int
    @test fin == 4

    #The vector form maps each frequency independently: 21 → 20 Hz (index 3),
    #41.5 → 40 Hz (index 5), 59 → 60 Hz (index 7)
    fv = [21.0, 41.5, 59.0]
    fvin = frqindex(fv,frqs)
    @test typeof(fvin) == Vector{Int}
    @test length(fvin) == length(fv)
    @test fvin[1] == 3
    @test fvin[2] == 5
    @test fvin[3] == 7
end
