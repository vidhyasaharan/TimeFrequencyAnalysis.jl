#Helper functions that read a number off a signal or a spectrum. They are tested against signals
#whose answer is known by construction, and - just as importantly - against the cases where the
#answer is not the one a caller might have wanted, since that is what the docstrings warn about.

@testset "peak_frequency" begin
    pfs = 1000.0
    n = 4000
    t = (0:n-1)./pfs

    #The dominant component of a tone in noise, from an array, a signal, or a spectrum already
    #computed - all three must agree
    x = 2.0.*cos.(2π*137 .*t) .+ 0.2.*randn(n)
    @test peak_frequency(x, pfs) ≈ 137 atol = 0.5
    @test peak_frequency(signal(x, pfs)) == peak_frequency(x, pfs)
    @test peak_frequency(magspec(signal(x, pfs))) == peak_frequency(x, pfs)

    #It finds the LARGEST component, not the lowest or the first: a strong low-frequency
    #interferer wins, which is the behaviour the docstring warns about
    y = cos.(2π*300 .*t) .+ 8.0.*cos.(2π*5 .*t)
    @test peak_frequency(y, pfs) ≈ 5 atol = 0.5

    #The answer is quantised to the DFT grid, so padding buys resolution. 137.4 Hz is off-grid at
    #ndft = 4000 (0.25 Hz bins) and lands closer with a finer grid.
    z = cos.(2π*137.4 .*t)
    coarse = peak_frequency(z, pfs; ndft = 1000)
    fine   = peak_frequency(z, pfs; ndft = 16000)
    @test abs(fine - 137.4) ≤ abs(coarse - 137.4)
    @test abs(fine - 137.4) < 0.1

    #Precision follows the signal, and a rectangular window is accepted like anywhere else
    @test peak_frequency(Float32.(x), 1000f0) isa Float32
    @test peak_frequency(x, pfs; wtype = "rect") ≈ 137 atol = 0.5
end
