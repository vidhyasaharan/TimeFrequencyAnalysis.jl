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
    t, h = gammatone_impulse_response(gf)

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
    t2, h2 = gammatone_impulse_response(gf; dur = 0.02)
    @test length(h2) == 320
end

@testset "filtering kernel" begin
    afs = 16000.0
    gf = gammatone_filter(afs, 1000.0)

    #A unit impulse through the kernel reproduces the exact binomial-form impulse response
    #(the closed form of gammatone_impulse_response), the strongest check the kernel can face
    N = 400
    e = zeros(N); e[1] = 1.0
    y = filt(gf, e)
    @test y isa Vector{ComplexF64}
    @test y[1] == gf.norm
    t, h = gammatone_impulse_response(gf; dur = N/afs)
    @test y ≈ h rtol = 1e-10

    #The measured forms agree with it too: filter_impresp through the kernel, and through
    #the expanded rational coefficients via DSP.filt — a third, independent filtering path
    @test filter_impresp(gf, N) ≈ h rtol = 1e-10
    @test filter_impresp(filter_coefs(gf), N) ≈ h rtol = 1e-8
    @test_throws ErrorException filter_impresp(gf, 0)

    #Linearity
    x1 = randn(500); x2 = randn(500)
    @test filt(gf, 2.0.*x1 .- 3.0.*x2) ≈ 2.0.*filt(gf, x1) .- 3.0.*filt(gf, x2)

    #Steady state: a unit tone at fc settles to envelope 1 (peak gain 2, halved because a
    #real tone splits across ±fc); off centre the envelope is |H(f)|/2, cross-validating
    #the kernel against the response machinery. Tolerance is loose: the filter's residual
    #response at the mirror frequency ripples the envelope at the 1e-3 level
    n = 0:7999
    for (f0, expected) in ((1000.0, 1.0), (1500.0, filter_magresp(gf, [1500.0])[1]/2))
        x = cos.(2π*f0.*n./afs)
        env = abs.(filt(gf, x))
        @test all(isapprox.(env[6001:8000], expected; rtol = 5e-3))
    end

    #A non-default order runs through the general (stage-by-stage) path and must match its
    #own closed form just as exactly
    gf6 = gammatone_filter(afs, 1000.0; order = 6)
    y6 = filt(gf6, e)
    _, h6 = gammatone_impulse_response(gf6; dur = N/afs)
    @test y6 ≈ h6 rtol = 1e-10

    #The in-place form matches, requires matching lengths, and allocates nothing after
    #warmup on both the order-4 and the general path
    x = randn(1000)
    for g in (gf, gf6)
        yb = Vector{ComplexF64}(undef, 1000)
        filt!(yb, g, x)
        @test yb ≈ filt(g, x)
        @test (@allocated filt!(yb, g, x)) == 0
    end
    @test_throws ErrorException filt!(Vector{ComplexF64}(undef, 5), gf, x)

    #The signal's precision decides the computation: a Float32 signal gives ComplexF32
    #output agreeing with the Float64 path, and a Float32-stored filter applied to a
    #Float64 signal stays Float64
    x32 = randn(Float32, 1000)
    y32 = filt(gf, x32)
    @test y32 isa Vector{ComplexF32}
    @test isapprox(y32, filt(gf, convert(Vector{Float64}, x32)); rtol = 1e-4)
    @test filt(gammatone_filter(16000f0, 1000f0), randn(100)) isa Vector{ComplexF64}
end

@testset "filterbank" begin
    afs = 16000.0
    fb = gammatone_filterbank(afs)
    @test fb isa gammatone_filterbank{Float64}
    @test length(fb.filters) == 30 == length(fb.fcs) #the gfb demo grid
    @test fb.fcs == erbfreq_array()
    @test 1000.0 ∈ fb.fcs
    @test all(f.fc == fc for (f,fc) in zip(fb.filters,fb.fcs))
    @test all(f.bw ≈ erb(f.fc) for f in fb.filters)
    @test all(f.order == 4 for f in fb.filters)

    #fmax defaults to min(6700, 0.45fs), so low sampling rates design cleanly; explicit
    #fmax at or above Nyquist fails with a clear error
    fb8 = gammatone_filterbank(8000.0)
    @test fb8.fcs[end] ≤ 0.45*8000
    @test_throws ErrorException gammatone_filterbank(afs; fmax = 9000.0)
    @test_throws ErrorException gammatone_filterbank(afs; fmin = 0.0)

    #Explicit centre frequencies, and the bank bw semantics: nothing → per-channel ERB
    #(above); scalar → the same bandwidth everywhere; vector → one per channel
    fcs = [200.0, 500.0, 1200.0]
    fbe = gammatone_filterbank(afs, fcs)
    @test fbe.fcs == fcs
    @test all(f.bw ≈ erb(f.fc) for f in fbe.filters)
    fbs = gammatone_filterbank(afs, fcs; bw = 100.0)
    @test all(f.bw == 100.0 for f in fbs.filters)
    fbv = gammatone_filterbank(afs, fcs; bw = [100.0, 150.0, 200.0])
    @test [f.bw for f in fbv.filters] == [100.0, 150.0, 200.0]
    @test_throws ErrorException gammatone_filterbank(afs, fcs; bw = [100.0, 150.0])
    @test all(f.bw ≈ 2*erb(f.fc) for f in gammatone_filterbank(afs, fcs; bw_scale = 2).filters)
    @test all(f.order == 3 for f in gammatone_filterbank(afs, fcs; order = 3).filters)
    @test_throws ErrorException gammatone_filterbank(afs, Float64[])

    @test gammatone_filterbank(16000f0) isa gammatone_filterbank{Float32}

    #Bank filtering: one row per channel, each row identical to the single-filter output;
    #the in-place form matches and checks its size
    x = randn(500)
    Y = filt(fb, x)
    @test Y isa Matrix{ComplexF64}
    @test size(Y) == (30, 500)
    @test Y[7,:] == filt(fb.filters[7], x)
    Yb = Matrix{ComplexF64}(undef, 30, 500)
    @test filt!(Yb, fb, x) == Y
    @test_throws ErrorException filt!(Matrix{ComplexF64}(undef, 29, 500), fb, x)

    #Measured bank impulse responses: one channel per row, identical to filtering an
    #impulse explicitly
    @test filter_impresp(fb, 100) == filt(fb, [1.0; zeros(99)])

    #Bank frequency responses: one row per channel, matching the single-filter methods
    fgrid = collect(500.0:50.0:2000.0)
    R = filter_resp(fb, fgrid)
    @test size(R) == (30, length(fgrid))
    @test R[7,:] == filter_resp(fb.filters[7], fgrid)
    M = filter_magresp(fb, fgrid)
    @test M[7,:] == filter_magresp(fb.filters[7], fgrid)
    @test M ≈ abs.(R)
end

@testset "analyses" begin
    #Containers on the shared fixture signal (fs = 8 kHz two-tone-in-noise from setup.jl)
    fb = gammatone_filterbank(fs)
    Y = gammatone_analysis(comp(), sig, fb)
    @test Y isa Matrix{ComplexF64}
    @test size(Y) == (length(fb.fcs), length(sig.x))

    tf = gammatone_analysis(sig, fb)
    @test tf isa timefreq{Float64, ComplexF64}
    @test tf.components == Y
    @test tf.frqs == fb.fcs
    @test tf.time == collect((0:length(sig.x)-1)./fs)
    @test tf.title == "Gammatone Filterbank (Hohmann 2002)"
    @test tf.frames === nothing

    cg = gammatone_cochleagram(sig, fb)
    @test cg isa timefreq{Float64, Float64}
    @test cg.components == abs.(Y) .+ eps(Float64)
    @test all(cg.components .> 0)
    @test cg.title == "Gammatone Cochleagram"
    @test gammatone_cochleagram(comp(), sig, fb) == cg.components

    #Keyword forms design the bank internally; the array form wraps the signal form
    tf2 = gammatone_analysis(sig)
    @test tf2.frqs == erbfreq_array(;fmax = 0.45*fs)
    @test gammatone_analysis(sig.x, fs).components == tf2.components
    @test gammatone_cochleagram(comp(), sig.x, fs) == gammatone_cochleagram(comp(), sig)

    #A tone at fcs[k] maximises the steady-state envelope in channel k
    fbt = gammatone_filterbank(16000.0)
    n = 0:3999
    for k in (5, 15, 25)
        xk = cos.(2π*fbt.fcs[k].*n./16000.0)
        E = abs.(filt(fbt, xk))
        @test argmax(vec(sum(E[:,2001:4000]; dims = 2))) == k
    end

    #Precision follows the signal, not the bank: a Float32 signal through the Float64 bank
    #gives Float32 containers throughout
    s32 = signal(randn(Float32, 2000), 8000f0)
    tf32 = gammatone_analysis(s32, fb)
    @test tf32 isa timefreq{Float32, ComplexF32}
    @test eltype(tf32.frqs) == Float32
    @test gammatone_cochleagram(s32, fb) isa timefreq{Float32, Float32}

    #A bank built for one sampling rate refuses a signal at another
    @test_throws ErrorException gammatone_analysis(signal(randn(100), 44100.0), fb)
end

@testset "delay helpers and summed response" begin
    afs = 16000.0
    gf = gammatone_filter(afs, 1000.0)
    γ = gf.order
    λ = Float64(abs(gf.coef))

    #Group delay at fc: the scalar form is the closed form γλ/(1-λ) samples; the numeric
    #curve agrees with it on a dense grid; and both sit within a few percent of the
    #analogue γ/(2πb) seconds (they differ at second order in b/fs)
    τfc = group_delay(gf)
    @test τfc ≈ γ*λ/(1 - λ)
    f = collect(800.0:1.0:1200.0)
    τ = group_delay(gf, f)
    @test length(τ) == length(f)
    @test τ[argmin(abs.(f .- 1000.0))] ≈ τfc rtol = 1e-3
    @test isapprox(τfc/afs, γ/(2π*Float64(gf.b)); rtol = 0.05)
    @test_throws ErrorException group_delay(gf, [1000.0])

    #The envelope peak sits at the exact discrete position ⌈(γλ-1)/(1-λ)⌉, a couple of
    #samples before the analogue prediction (γ-1)/(2πb) — the discrete binomial envelope
    #peaks earlier than its continuous counterpart
    for gfx in (gf, gammatone_filter(afs, 100.0), gammatone_filter(afs, 1000.0, 200.0, 3.0))
        γx = gfx.order
        λx = Float64(abs(gfx.coef))
        nd = envelope_delay(gfx)
        @test nd == ceil(Int, (γx*λx - 1)/(1 - λx))
        @test abs(nd - (γx - 1)*afs/(2π*Float64(gfx.b))) ≤ 4
    end

    #filter_resp equals the DTFT of the (fully decayed) impulse response: an independent
    #time-domain path to the same frequency response
    _, h = gammatone_impulse_response(gf; dur = 0.25)
    for f0 in (250.0, 1000.0, 2000.0)
        θ0 = 2π*f0/afs
        Hdtft = sum(h[n+1]*cis(-θ0*n) for n ∈ 0:length(h)-1)
        @test Hdtft ≈ filter_resp(gf, [f0])[1] rtol = 1e-6
    end

    #Summed bank response: the channel-wise sum, reducing to the single response for a
    #one-channel bank, and far from flat without alignment (the step-6 motivation)
    fb = gammatone_filterbank(afs)
    fgrid = collect(100.0:5.0:7900.0)
    sr = summed_resp(fb, fgrid)
    @test sr ≈ sum(filter_resp(gfx, fgrid) for gfx in fb.filters)
    fb1 = gammatone_filterbank(afs, [1000.0])
    @test summed_resp(fb1, fgrid) ≈ filter_resp(fb1.filters[1], fgrid)
    inband = abs.(sr[frqindex(300.0, fgrid):frqindex(5000.0, fgrid)])
    @test maximum(inband)/minimum(inband) > 1.5
end

@testset "delay and phase compensation" begin
    afs = 16000.0

    #Cross-check against the reference recipe (gfb_delay_new / hohmann2002delay, as ported
    #by pyfar and pyfilterbank) run on the actual pyfilterbank filters: a 3-channel
    #Route-B bank (bw = erb(fc), edges -3 dB) aligned to 4 ms. Reference values generated
    #on 2026-08-14 by "scripts/generate-gammatone-refs.py". The 250 Hz channel is too slow
    #for the window and clamps at its edge; the other two align.
    fbB = gammatone_filterbank(afs, [250.0, 1000.0, 3000.0]; attenuation_db = 3)
    dB = gammatone_delay(fbB; delay = 0.004)
    @test dB.nd == 64
    @test dB.delays == [0, 16, 47]
    for (k, ref) in enumerate((0.9757621274052979 + 0.21883388842107543im,
                               0.9999999920136321 + 0.00012638328942236236im,
                               0.3830361826798044 - 0.9237333396376268im))
        @test isapprox(dB.phase_factors[k], ref; rtol = 1e-9)
    end
    @test all(abs.(dB.phase_factors) .≈ 1)

    #Alignment on the default bank: after compensation every channel that peaks inside
    #the window has its envelope maximum exactly at nd and is real-positive there, the
    #summed real part pulses at nd, and the slow bottom channels clamp to Δ = 0
    fb = gammatone_filterbank(afs)
    d = gammatone_delay(fb) #default 4 ms → nd = 64
    @test d.nd == 64
    e = zeros(300); e[1] = 1.0
    Y = filt(fb, e)
    Yc = compensate(Y, d)
    @test size(Yc) == size(Y)
    for k in eachindex(fb.filters)
        if envelope_delay(fb.filters[k]) ≤ d.nd
            row = view(Yc, k, :)
            @test argmax(abs.(row)) == d.nd + 1
            @test real(row[d.nd + 1]) > 0.999*abs(row[d.nd + 1])
        else
            @test d.delays[k] == 0
        end
    end
    @test any(d.delays .== 0)
    @test argmax(vec(sum(real.(Yc); dims = 1))) == d.nd + 1

    #compensate copies, compensate! matches it in place, and shapes/arguments are checked
    @test Yc != Y
    Yb = copy(Y)
    @test compensate!(Yb, d) == Yc
    @test_throws ErrorException compensate(Y[1:5, :], d)
    @test_throws ErrorException gammatone_delay(fb; delay = -0.001)

    #timefreq form preserves the axes and marks the title; the align keyword is the same
    #compensation folded into the analyses; the cochleagram container refuses (no phase)
    s = signal(randn(600), afs)
    tf = gammatone_analysis(s, fb)
    tfc = compensate(tf, d)
    @test tfc isa timefreq{Float64, ComplexF64}
    @test tfc.components == compensate(tf.components, d)
    @test tfc.frqs == tf.frqs && tfc.time == tf.time
    @test endswith(tfc.title, "(delay compensated)")
    @test_throws ErrorException compensate(gammatone_cochleagram(s, fb), d)
    @test gammatone_analysis(s, fb; align = 0.004).components == tfc.components
    @test gammatone_analysis(s, fb; align = 0.004).title == tfc.title
    @test gammatone_cochleagram(comp(), s, fb; align = 0.004) == abs.(tfc.components) .+ eps(Float64)
    @test endswith(gammatone_cochleagram(s, fb; align = 0.004).title, "(delay compensated)")
    @test gammatone_analysis(s; align = 0.004) isa timefreq{Float64, ComplexF64}

    #Precision: a Float64-built compensation applies to Float32 subbands in Float32
    s32 = signal(randn(Float32, 400), 16000f0)
    tfc32 = compensate(gammatone_analysis(s32, fb), d)
    @test tfc32 isa timefreq{Float32, ComplexF32}

    #Compensated summed response: the closed-form sum Σ c_k e^{-iθΔ_k} H_k, accepting a
    #prebuilt object or a target delay. What alignment delivers is PHASE coherence: after
    #removing the pure nd-sample delay, the compensated sum's phase is nearly flat where
    #the channels are aligned, while the raw sum's phase swings the full ±π. (Magnitude
    #flatness is NOT improved by alignment — the incoherent raw sum is statistically
    #smooth while the coherent sum shows the picket-fence ripple the v2 mixer removes —
    #so only a loose magnitude bound is asserted.)
    fgrid = collect(100.0:5.0:7900.0)
    srr = summed_resp(fb, fgrid)
    src = summed_resp(fb, fgrid; delay = d)
    @test src ≈ summed_resp(fb, fgrid; delay = 0.004)
    θ = 2π.*fgrid./afs
    manual = sum(d.phase_factors[k].*cis.(-θ.*d.delays[k]).*filter_resp(fb.filters[k], fgrid) for k in eachindex(fb.filters))
    @test src ≈ manual
    aligned_band = frqindex(800.0, fgrid):frqindex(5000.0, fgrid)
    @test maximum(abs.(angle.(src[aligned_band].*cis.(θ[aligned_band].*d.nd)))) < 0.2
    @test maximum(abs.(angle.(srr[aligned_band].*cis.(θ[aligned_band].*d.nd)))) > 2
    band = frqindex(300.0, fgrid):frqindex(5000.0, fgrid)
    @test maximum(abs.(src[band]))/minimum(abs.(src[band])) < 2.5
end

#default_align exists because a fixed alignment target cannot serve every bank: a channel too
#slow to peak inside the search window clamps to Delta = 0 and stays misaligned. The testset
#above asserts the default 4 ms bank does exactly that (any(d.delays .== 0)); here we assert
#:auto does not. compensation_lead/trim then remove the samples compensation invents.
@testset "automatic alignment target and lead trimming" begin
    afs = 16000.0
    fb = gammatone_filterbank(afs)

    #The bank method is the per-filter one applied across channels, and delays fall with fc
    @test envelope_delay(fb) == [envelope_delay(gf) for gf in fb.filters]
    @test envelope_delay(fb) isa Vector{Int}
    @test issorted(envelope_delay(fb); rev = true)

    #The target is one sample past the slowest channel, so nd exceeds every envelope delay
    ta = default_align(fb)
    @test ta isa Float64
    @test ta ≈ (maximum(envelope_delay(fb)) + 1)/afs
    da = gammatone_delay(fb; delay = :auto)
    @test da.nd == maximum(envelope_delay(fb)) + 1
    @test da.nd == gammatone_delay(fb; delay = ta).nd

    #The point of :auto — every channel peaks strictly inside the window, so none is clamped.
    #The 4 ms default on this same bank does clamp, which is the bug being fixed.
    @test all(envelope_delay(fb) .< da.nd + 1)
    @test all(da.delays .> 0)
    @test any(gammatone_delay(fb; delay = 0.004).delays .== 0)

    #With :auto every channel's impulse response peaks at the common target instant
    e = zeros(1200); e[1] = 1.0
    Ya = compensate(filt(fb, e), da)
    for k in eachindex(fb.filters)
        @test argmax(abs.(view(Ya, k, :))) == da.nd + 1
    end

    #compensation_lead is the number of leading columns holding zero-fill in some channel
    @test compensation_lead(da) == maximum(da.delays)
    @test all(Ya[argmax(da.delays), 1:compensation_lead(da)] .== 0)

    #trim drops exactly those columns from every channel, keeping the matrix rectangular and
    #the channels aligned with one another
    Yt = compensate(filt(fb, e), da; trim = true)
    @test size(Yt) == (size(Ya, 1), size(Ya, 2) - compensation_lead(da))
    @test Yt == Ya[:, (compensation_lead(da) + 1):end]
    @test compensate(filt(fb, e), da; trim = false) == Ya
    #no channel now starts on an invented zero
    @test all(abs.(Yt[:, 1]) .> 0)

    #timefreq form trims the time axis to match and keeps the original origin rather than
    #re-zeroing it, so the result still refers to the same instants as the input signal
    s = signal(randn(1200), afs)
    tf = gammatone_analysis(s, fb)
    tft = compensate(tf, da; trim = true)
    @test size(tft.components, 2) == length(tft.time)
    @test tft.time == tf.time[(compensation_lead(da) + 1):end]
    @test tft.time[1] ≈ compensation_lead(da)/afs
    @test tft.time[1] > 0
    @test compensate(tf, da).time == tf.time      #untrimmed still starts at 0
    @test endswith(tft.title, "(delay compensated)")

    #Errors: a bad symbol, and a signal shorter than the lead being trimmed
    @test_throws ErrorException gammatone_delay(fb; delay = :nearest)
    @test_throws ErrorException compensate(filt(fb, zeros(10)), da; trim = true)
    @test_throws ErrorException compensate(gammatone_cochleagram(s, fb), da; trim = true)

    #Precision follows the bank, and :auto composes with the align keyword on the analyses
    fb32 = gammatone_filterbank(16000f0)
    @test default_align(fb32) isa Float32
    @test gammatone_analysis(s, fb; align = ta).components == compensate(tf, da).components
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
    t32, h32 = gammatone_impulse_response(gf32)
    @test eltype(t32) == Float32
    @test eltype(h32) == ComplexF32
end
