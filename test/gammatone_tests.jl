#Gammatone filter design (Hohmann 2002), checked three ways: closed-form invariants (peak
#gain exactly 2 at fc, a₄ = 5π/16, the binomial impulse response), behavioural checks of
#the bandwidth semantics (ERB default, explicit bw, bw_scale, edge attenuation), and
#cross-implementation reference values. The reference literals were generated on 2026-08-13
#by "scripts/generate-gammatone-refs.py": Route B values from the actual pyfilterbank
#implementation (pyfilterbank/gammatone.py, github.com/SiggiGue/pyfilterbank master,
#design_filter with attenuation_half_bandwidth_db = -3, Python 3.9/NumPy 2.0); Route A
#values from an independent line-by-line transcription of Hohmann (2002) Eqs. 13-17
#(pyfilterbank does not implement that route - its a_gamma branch is commented out - and
#AMT's hohmann2002filter implements the same equations).

@testset "design (Route A: ERB matching)" begin
    afs = 16000.0
    gf = gammatone_filter(afs, 1000.0)
    @test gf isa gammatone_filter{Float64}
    @test gf.order == 4
    @test 0 < abs(gf.coef) < 1
    @test gf.bw ≈ erb(1000.0)
    @test gf.b ≈ gf.bw/(5π/16) #a₄ = 5π/16 exactly
    @test TimeFrequencyAnalysis.a_gamma(4) ≈ 5π/16

    #The magnitude response peaks exactly at fc with gain exactly 2, and the response at fc
    #is real and positive (the cascade denominator (1-λ)^order is real there)
    @test filter_resp(gf, [1000.0])[1] ≈ 2
    fgrid = collect(0.0:1.0:8000.0)
    mres = filter_magresp(gf, fgrid)
    @test fgrid[argmax(mres)] == 1000.0

    #The expanded-polynomial evaluation path must agree with the closed form of the cascade
    for f in (500.0, 1000.0, 2500.0)
        Hcf = gf.norm/(1 - gf.coef*exp(-im*2π*f/afs))^gf.order
        @test filter_resp(gf, [f])[1] ≈ Hcf
    end

    #Cross-check against Hohmann (2002) Eqs. 13-17 reference values (provenance above)
    gf100 = gammatone_filter(afs, 100.0)
    gf4k = gammatone_filter(afs, 4000.0)
    gfhalf = gammatone_filter(afs, 1000.0; bw_scale = 0.5)
    gf3rd = gammatone_filter(44100.0, 2000.0; bw_scale = 2, order = 3)
    for (gfx, coef, norm, b) in (
            (gf,     0.8761422420079817 + 0.3629099992076765im,  1.4256065230178252e-5, 135.09894743398627),
            (gf100,  0.9851428880070562 + 0.03870636948764195im, 7.898394665881612e-8,  36.153186806968755),
            (gf4k,   5.1014124748941916e-17 + 0.8331238816687381im, 0.001550981981190936, 464.91814952404457),
            (gfhalf, 0.899694328624853 + 0.3726655929065704im,   9.392206309200397e-7,  67.54947371699313),
            (gf3rd,  0.9054283035531185 + 0.2652209322527693im,  0.0003612289003736513, 408.3978024400095))
        @test isapprox(gfx.coef, coef; rtol = 1e-12)
        @test isapprox(gfx.norm, norm; rtol = 1e-12)
        @test isapprox(gfx.b, b; rtol = 1e-12)
    end
end

@testset "bandwidth semantics" begin
    afs = 16000.0

    #Explicit bw with the ERB value reproduces the default exactly; bw_scale composes
    #multiplicatively with both the default and an explicit bw
    gf = gammatone_filter(afs, 1000.0)
    @test gammatone_filter(afs, 1000.0; bw = erb(1000.0)).coef == gf.coef
    @test gammatone_filter(afs, 1000.0; bw_scale = 2).coef == gammatone_filter(afs, 1000.0; bw = 2*erb(1000.0)).coef
    gfc = gammatone_filter(afs, 1000.0; bw = 200.0, bw_scale = 1.5)
    @test gfc.bw ≈ 300.0
    @test gfc.coef == gammatone_filter(afs, 1000.0; bw = 300.0).coef

    #Wider bandwidth means faster envelope decay: |ã| falls monotonically with bw
    λs = [abs(gammatone_filter(afs, 1000.0; bw = bwv).coef) for bwv in (50.0, 132.6, 300.0)]
    @test issorted(λs; rev = true)
end

@testset "design (Route B: edge attenuation)" begin
    afs = 16000.0

    #Band edges fc ± bw/2 sit exactly attenuation_db below the peak
    for A in (3.0, 6.0)
        gfb = gammatone_filter(afs, 1000.0, 200.0, A)
        medge = filter_magresp(gfb, [900.0, 1100.0])
        mpeak = filter_magresp(gfb, [1000.0])[1]
        @test medge[1]/mpeak ≈ 10.0^(-A/20) rtol = 1e-6
        @test medge[2]/mpeak ≈ 10.0^(-A/20) rtol = 1e-6
    end

    #Cross-check against the actual pyfilterbank design_filter (provenance above)
    for (gfx, coef, norm) in (
            (gammatone_filter(16000.0, 1000.0, 200.0, 3.0),
                0.8440145543437426 + 0.34960227524946186im, 0.00011168481292224301),
            (gammatone_filter(16000.0, 1000.0, erb(1000.0), 3.0), #pyfilterbank's band_width_factor=1 path
                0.8701016253738603 + 0.36040789387272687im, 2.296065111091373e-5),
            (gammatone_filter(44100.0, 4000.0, 456.5, 3.0),
                0.7812092339006116 + 0.500624121205301im, 5.4186074022075584e-5),
            (gammatone_filter(16000.0, 500.0, 100.0, 3.0; order = 2),
                0.9512577418163434 + 0.18921692941291635im, 0.001812744551343129))
        @test isapprox(gfx.coef, coef; rtol = 1e-12)
        @test isapprox(gfx.norm, norm; rtol = 1e-12)
    end

    #For the same nominal bandwidth, pinning the -3 dB width (Route B) gives a wider filter
    #(smaller pole radius) than matching the equivalent rectangular bandwidth (Route A)
    @test abs(gammatone_filter(afs, 1000.0, erb(1000.0), 3.0).coef) < abs(gammatone_filter(afs, 1000.0).coef)
end

@testset "filter_coefs conversion" begin
    gf = gammatone_filter(16000.0, 1000.0)
    F = filter_coefs(gf)
    @test F isa filter_coefs{Float64, ComplexF64}
    @test F.num == [gf.norm]
    @test length(F.den) == 5
    @test F.den[1] == 1
    @test F.den[2] ≈ -4*gf.coef
    @test F.den[5] ≈ gf.coef^4
end

@testset "impulse response (closed form)" begin
    afs = 16000.0
    gf = gammatone_filter(afs, 1000.0)
    t, h = impulse_response(gf)

    #Time axis and the first terms of h[n] = k⋅binomial(n+3,3)⋅ãⁿ
    @test t[1] == 0
    @test t[2] ≈ 1/afs
    @test length(t) == length(h)
    @test h[1] ≈ gf.norm
    @test h[2] ≈ 4*gf.norm*gf.coef
    @test h[11] ≈ gf.norm*binomial(13,3)*gf.coef^10

    #Default duration runs until the envelope has decayed to -80 dB re its peak; the
    #envelope peak sits near the analogue prediction (γ-1)/(2πb)
    env = abs.(h)
    @test env[end] ≤ 1e-4*maximum(env)*(1 + 1e-10)
    npeak_analog = (gf.order - 1)*afs/(2π*gf.b)
    @test abs((argmax(env) - 1) - npeak_analog) ≤ 3

    #dur keyword overrides the length
    t2, h2 = impulse_response(gf; dur = 0.02)
    @test length(h2) == 320
end

@testset "filtering kernel" begin
    afs = 16000.0
    gf = gammatone_filter(afs, 1000.0)

    #A unit impulse through the kernel reproduces the exact binomial-form impulse response
    #(the closed form of impulse_response), the strongest check the kernel can face
    N = 400
    e = zeros(N); e[1] = 1.0
    y = gammatone_filt(gf, e)
    @test y isa Vector{ComplexF64}
    @test y[1] == gf.norm
    t, h = impulse_response(gf; dur = N/afs)
    @test y ≈ h rtol = 1e-10

    #Linearity
    x1 = randn(500); x2 = randn(500)
    @test gammatone_filt(gf, 2.0.*x1 .- 3.0.*x2) ≈ 2.0.*gammatone_filt(gf, x1) .- 3.0.*gammatone_filt(gf, x2)

    #Steady state: a unit tone at fc settles to envelope 1 (peak gain 2, halved because a
    #real tone splits across ±fc); off centre the envelope is |H(f)|/2, cross-validating
    #the kernel against the response machinery. Tolerance is loose: the filter's residual
    #response at the mirror frequency ripples the envelope at the 1e-3 level
    n = 0:7999
    for (f0, expected) in ((1000.0, 1.0), (1500.0, filter_magresp(gf, [1500.0])[1]/2))
        x = cos.(2π*f0.*n./afs)
        env = abs.(gammatone_filt(gf, x))
        @test all(isapprox.(env[6001:8000], expected; rtol = 5e-3))
    end

    #A non-default order runs through the general (stage-by-stage) path and must match its
    #own closed form just as exactly
    gf6 = gammatone_filter(afs, 1000.0; order = 6)
    y6 = gammatone_filt(gf6, e)
    _, h6 = impulse_response(gf6; dur = N/afs)
    @test y6 ≈ h6 rtol = 1e-10

    #The in-place form matches, requires matching lengths, and allocates nothing after
    #warmup on both the order-4 and the general path
    x = randn(1000)
    for g in (gf, gf6)
        yb = Vector{ComplexF64}(undef, 1000)
        gammatone_filt!(yb, g, x)
        @test yb ≈ gammatone_filt(g, x)
        @test (@allocated gammatone_filt!(yb, g, x)) == 0
    end
    @test_throws ErrorException gammatone_filt!(Vector{ComplexF64}(undef, 5), gf, x)

    #The signal's precision decides the computation: a Float32 signal gives ComplexF32
    #output agreeing with the Float64 path, and a Float32-stored filter applied to a
    #Float64 signal stays Float64
    x32 = randn(Float32, 1000)
    y32 = gammatone_filt(gf, x32)
    @test y32 isa Vector{ComplexF32}
    @test isapprox(y32, gammatone_filt(gf, convert(Vector{Float64}, x32)); rtol = 1e-4)
    @test gammatone_filt(gammatone_filter(16000f0, 1000f0), randn(100)) isa Vector{ComplexF64}
end

@testset "argument checking and precision" begin
    afs = 16000.0
    @test_throws ErrorException gammatone_filter(afs, 9000.0)  #fc above fs/2
    @test_throws ErrorException gammatone_filter(afs, 0.0)     #fc must be positive
    @test_throws ErrorException gammatone_filter(afs, 1000.0; bw = -10.0)
    @test_throws ErrorException gammatone_filter(afs, 1000.0; order = 11)
    @test_throws ErrorException gammatone_filter(afs, 1000.0, 200.0, -3.0) #attenuation must be positive
    @test_logs (:warn, r"Nyquist") gammatone_filter(afs, 7500.0)

    #Integer inputs promote to Float64; Float32 inputs stay Float32 end to end, matching
    #the Float64 design up to Float32 rounding
    @test gammatone_filter(16000, 1000) isa gammatone_filter{Float64}
    gf32 = gammatone_filter(16000f0, 1000f0)
    @test gf32 isa gammatone_filter{Float32}
    @test gf32.coef isa ComplexF32
    @test abs(gf32.coef - gammatone_filter(16000.0, 1000.0).coef) < 1e-6
    t32, h32 = impulse_response(gf32)
    @test eltype(t32) == Float32
    @test eltype(h32) == ComplexF32
end
