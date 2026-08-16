#Digital filter representation and frequency responses, checked against filters whose
#responses are known in closed form (identity, two-tap moving average, single-pole
#lowpass, one-pole complex resonator), plus the frequency-conversion helpers and
#coefficient type promotion (real and complex).

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

@testset "complex coefficients" begin
    afs = 8000.0

    #Promotion rules: any complex input promotes both vectors to Complex of the common real
    #precision; the first type parameter remains the real precision, so the {T} supertype
    #checks used for real filters keep working
    Fc = filter_coefs([1.0], [1.0, 0.5im])
    @test Fc isa filter_coefs{Float64, ComplexF64}
    @test Fc isa filter_coefs{Float64}
    @test filter_coefs(ComplexF32[1], ComplexF32[1, 0.5]) isa filter_coefs{Float32, ComplexF32}
    @test filter_coefs(Float32[1], ComplexF32[1, 0.5]) isa filter_coefs{Float32, ComplexF32} #real/complex mix
    @test filter_coefs(Float32[1], ComplexF64[1, 0.5]) isa filter_coefs{Float64, ComplexF64} #mixed precision promotes
    @test filter_coefs([1, 2], [1, -im]) isa filter_coefs{Float64, ComplexF64} #integer input promotes

    #One-pole complex resonator 1/(1 - ã z⁻¹) with ã = λ e^{iβ}: peak gain 1/(1-λ) exactly
    #at the pole frequency, and no mirror peak at -fp — the response is one-sided, which no
    #real-coefficient filter can achieve
    λ = 0.9
    fp = 1000.0
    ã = λ*exp(im*2π*fp/afs)
    Fp = filter_coefs([1.0], [1.0, -ã])
    @test filter_magresp(Fp, [fp], afs)[1] ≈ 1/(1 - λ)
    fgrid = collect(-4000.0:10.0:4000.0)
    mres = filter_magresp(Fp, fgrid, afs)
    @test fgrid[argmax(mres)] == fp
    @test filter_magresp(Fp, [-fp], afs)[1] < 1 #mirror frequency sits far below the peak gain of 10

    #The two evaluation paths must agree for complex coefficients too
    @test filter_magresp(Fp, fgrid, afs) ≈ abs.(filter_resp(Fp, fgrid, afs))
end

@testset "filt and impulse responses" begin
    #filt on a filter_coefs is DSP's coefficient filt; the in-place form matches
    Fap0 = filter_coefs([1], [1.0, -0.9])
    x = randn(200)
    @test filt(Fap0, x) == filt(Fap0.num, Fap0.den, x)
    yb = zeros(200)
    @test filt!(yb, Fap0, x) == filt(Fap0, x)

    #filter_impresp measures by running a unit impulse through any filter with a filt
    #method: FIR coefficients reproduce themselves, a real one-pole gives 0.9ⁿ, a complex
    #one-pole gives ãⁿ; like the frequency-response helpers the measurement runs in Float64
    Ffir = filter_coefs([1.0, 2.0, 3.0], [1.0])
    @test filter_impresp(Ffir, 5) == [1.0, 2.0, 3.0, 0.0, 0.0]
    Fap = filter_coefs([1], [1.0, -0.9])
    @test filter_impresp(Fap, 20) ≈ 0.9.^(0:19)
    ã = 0.9*exp(im*2π*1000.0/8000.0)
    Fp = filter_coefs([1.0], [1.0, -ã])
    h = filter_impresp(Fp, 20)
    @test h isa Vector{ComplexF64}
    @test h ≈ ã.^(0:19)
    @test_throws ErrorException filter_impresp(Fp, 0)
end

@testset "z⁻¹ phase convention" begin
    afs = 8000.0

    #A pure one-sample delay (num = [0,1]) has response e^{-iθ}: phase -2πf/fs, magnitude 1
    Fd = filter_coefs([0.0, 1.0], [1.0])
    for f in (500.0, 1000.0, 3000.0)
        @test angle(TimeFrequencyAnalysis.H(Fd, f, afs)) ≈ -2π*f/afs
        @test TimeFrequencyAnalysis.Hmag(Fd, f, afs) ≈ 1
    end

    #Real-coefficient filters keep conjugate-symmetric responses: H(-f) = conj(H(f))
    Fap = filter_coefs([1], [1.0, -0.9])
    @test TimeFrequencyAnalysis.H(Fap, -1000.0, afs) ≈ conj(TimeFrequencyAnalysis.H(Fap, 1000.0, afs))
end
