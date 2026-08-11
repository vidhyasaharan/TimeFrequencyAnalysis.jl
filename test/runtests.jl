#One test file per source file. setup.jl provides the shared fixtures (the fs/x/signal
#test signal) and seeds the RNG so that the stochastic tests are deterministic.
include("setup.jl")

@testset verbose = true "TimeFrequencyAnalysis" begin
    @testset "types" begin
        include("types_tests.jl")
    end
    @testset "sampling" begin
        include("sampling_tests.jl")
    end
    @testset "framing" begin
        include("framing_tests.jl")
    end
    @testset "windows" begin
        include("windows_tests.jl")
    end
    @testset "freq. grids" begin
        include("grids_tests.jl")
    end
    @testset "generators" begin
        include("generators_tests.jl")
    end
    @testset "spectral" begin
        include("spectral_tests.jl")
    end
    @testset "correlations" begin
        include("correlations_tests.jl")
    end
    @testset "filters" begin
        include("filters_tests.jl")
    end
    @testset "element ops" begin
        include("element_ops_tests.jl")
    end
    @testset "recipes" begin
        include("recipes_tests.jl")
    end
    @testset "Float32" begin
        include("float32_tests.jl")
    end
end
