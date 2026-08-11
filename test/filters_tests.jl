#Digital filter representation and frequency responses, checked against filters whose
#responses are known in closed form (identity, two-tap moving average, single-pole
#lowpass), plus the frequency-conversion helpers and coefficient type promotion.

@testset "frequency conversion" begin
    #freq2θ maps Hz to radians per sample (DC → 0, Nyquist → π); freq2z places those
    #frequencies on the unit circle (DC → 1, Nyquist → -1, always magnitude 1)
    afs = 8000.0
    @test TimeFrequencyAnalysis.freq2θ(0, afs) == 0
    @test TimeFrequencyAnalysis.freq2θ(afs/2, afs) ≈ pi
    f = [0.0, 1000.0, 4000.0]
    @test TimeFrequencyAnalysis.freq2θ(f, afs) ≈ 2pi .* f ./ afs
    z = TimeFrequencyAnalysis.freq2z(f, afs)
    @test all(abs.(z) .≈ 1) #frequencies map onto the unit circle
    @test z[1] ≈ 1
    @test TimeFrequencyAnalysis.freq2z(afs/2, afs) ≈ -1
end

@testset "filter_coefs" begin
    #The convenience constructor must promote both coefficient vectors to one common
    #floating point element type (integers to Float64, mixed Float32/Float64 to Float64)
    F = filter_coefs([1, 2], [1, 0, 1]) #integer input promotes to Float64
    @test F isa filter_coefs{Float64}
    @test F.num == [1.0, 2.0]
    @test F.den == [1.0, 0.0, 1.0]
    @test filter_coefs(Float32[1], Float32[1, -0.5]) isa filter_coefs{Float32}
    @test filter_coefs(Float32[1], [1.0, -0.5]) isa filter_coefs{Float64} #mixed promotes
end

@testset "responses" begin
    afs = 8000.0
    f = collect(0.0:100.0:4000.0)

    #powers is the polynomial-evaluation helper: [1, x, x², x³]
    @test TimeFrequencyAnalysis.powers(2.0, 3) == [1.0, 2.0, 4.0, 8.0]

    #Identity filter (num = den = [1]): unit response at every frequency
    Fid = filter_coefs([1],[1])
    @test filter_resp(Fid, f, afs) ≈ ones(length(f))
    @test filter_magresp(Fid, f, afs) ≈ ones(length(f))

    #Two-tap moving average [1,1]: |H| = |1 + e^(iθ)| is 2 at DC and 0 at Nyquist
    Ffir = filter_coefs([1,1],[1])
    m = filter_magresp(Ffir, [0.0, afs/2], afs)
    @test m[1] ≈ 2
    @test m[2] ≈ 0 atol = 1e-12

    #The two evaluation paths must agree: the optimised Hmag equals abs of the complex
    #response computed via H/powers, pointwise and over a grid
    Fap = filter_coefs([1], [1.0, -0.9])
    @test filter_magresp(Fap, f, afs) ≈ abs.(filter_resp(Fap, f, afs))
    @test TimeFrequencyAnalysis.H(Fap, 1000.0, afs) ≈ filter_resp(Fap, [1000.0], afs)[1]
    @test TimeFrequencyAnalysis.Hmag(Fap, 1000.0, afs) ≈ abs(TimeFrequencyAnalysis.H(Fap, 1000.0, afs))

    #Single-pole lowpass 1/(1 - 0.9z): magnitude must decrease monotonically from DC to
    #Nyquist over [0, fs/2]
    mres = filter_magresp(Fap, f, afs)
    @test all(diff(mres) .< 0)
end
