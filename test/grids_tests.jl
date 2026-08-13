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

@testset "ERB scale" begin
    #Glasberg & Moore constants in the Hohmann (2002) form: ERB(0) = 24.7 Hz and
    #ERB(1 kHz) ≈ 132.63 Hz; the ERB-rate scale is 0 at 0 Hz and inverts exactly
    @test erb(0) == 24.7
    @test erb(1000.0) ≈ 132.6331 atol = 1e-3
    @test freq2erb(0) == 0
    @test erb2freq(0) == 0
    for f in (50.0, 250.0, 1000.0, 4000.0)
        @test erb2freq(freq2erb(f)) ≈ f
    end
end

@testset "erbfreq_array" begin
    #The default grid is the ~30-channel gfb demo grid (70 Hz to 6.7 kHz around 1 kHz at
    #one filter per ERB): sorted, inside [fmin, fmax], uniformly spaced on the ERB-rate
    #scale, with the anchor present exactly
    g = erbfreq_array()
    @test length(g) == 30
    @test issorted(g)
    @test g[1] ≥ 70 && g[end] ≤ 6700
    @test 1000.0 ∈ g #pinned exactly, not just approximately
    E = freq2erb.(g)
    @test all(isapprox.(diff(E), 1.0; atol = 1e-8))

    #Doubling the density halves the ERB spacing; the anchor is pinned for any base_frq
    g2 = erbfreq_array(;filters_per_erb = 2)
    @test all(isapprox.(diff(freq2erb.(g2)), 0.5; atol = 1e-8))
    @test 432.1 ∈ erbfreq_array(;fmin = 100, fmax = 2000, base_frq = 432.1)

    #The anchor must lie inside [fmin, fmax] and the density must be positive
    @test_throws ErrorException erbfreq_array(;fmin = 100, fmax = 2000, base_frq = 50)
    @test_throws ErrorException erbfreq_array(;filters_per_erb = 0)
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
