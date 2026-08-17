#Modulation analysis. The invariants here are analytic: an amplitude modulated tone has an
#envelope whose mean, extremes and harmonic structure are all known in closed form, so the
#envelope stage can be checked against theory rather than against another implementation.
#
#fs = 3125 Hz throughout, matching the sonar case the defaults were chosen for, so setup.jl's
#8 kHz signal is not used.

const mfs = 3125.0

#A bank spanning a decade around 200 Hz, wide enough (bw_scale = 2) that a 7 Hz modulation's
#sidebands sit well inside every channel
mbank(; fmin = 100, fmax = 400, base_frq = 200, filters_per_erb = 2, kwargs...) =
    gammatone_filterbank(mfs; fmin, fmax, base_frq, filters_per_erb, bw_scale = 2, kwargs...)

am_tone(A, m, fc, fm, dur; fs = mfs) =
    (t = (0:round(Int, dur*fs)-1)./fs; A.*(1 .+ m.*cos.(2π*fm.*t)).*cos.(2π*fc.*t))

@testset "envelope amplitude and shape" begin
    fb = mbank()
    k = frqindex(200.0, fb.fcs)
    @test fb.fcs[k] ≈ 200

    #The gain-2 normalisation makes the envelope read in the units of the signal: an unmodulated
    #tone of amplitude A gives an amplitude envelope of exactly A and a power envelope of A^2.
    #This is what lets a modulation index be read off without calibration.
    for A in (1.0, 0.3)
        s = signal(A.*cos.(2π*200 .*(0:round(Int,20mfs)-1)./mfs), mfs)
        ea = power_envelope(comp(), s, fb; envelope = :amplitude)[k,:]
        ep = power_envelope(comp(), s, fb; envelope = :power)[k,:]
        @test all(isapprox.(ea, A; rtol = 1e-3))
        @test all(isapprox.(ep, A^2; rtol = 1e-3))
        #and nothing of the filter start-up survives the default trim
        @test maximum(abs.(ep .- A^2)) < 1e-3*A^2
    end

    #Amplitude modulated: the power envelope is A^2(1 + m cos 2*pi*fm*t)^2, so its mean is
    #A^2(1 + m^2/2) and it swings between A^2(1-m)^2 and A^2(1+m)^2
    A, m = 1.0, 0.5
    s = signal(am_tone(A, m, 200.0, 7.0, 20.0), mfs)
    e = power_envelope(comp(), s, fb; envelope = :power)[k,:]
    @test sum(e)/length(e) ≈ A^2*(1 + m^2/2) rtol = 0.01
    @test minimum(e) ≈ A^2*(1-m)^2 rtol = 0.05
    @test maximum(e) ≈ A^2*(1+m)^2 rtol = 0.05

    #The amplitude envelope of the same signal is A(1 + m cos), i.e. the square root of the above.
    #Undecimated the two are exactly related, both being taken from the same complex subband...
    eu_a = power_envelope(comp(), s, fb; envelope = :amplitude, fs_env = mfs)[k,:]
    eu_p = power_envelope(comp(), s, fb; envelope = :power, fs_env = mfs)[k,:]
    @test eu_a ≈ sqrt.(eu_p) rtol = 1e-12

    #...but decimation is a LINEAR operation and so does not commute with squaring: the decimated
    #amplitude envelope is only approximately the square root of the decimated power envelope.
    #Choosing :power or :amplitude therefore has to happen before the envelope is decimated, which
    #is why it is a keyword on power_envelope rather than something a caller can apply afterwards.
    ea = power_envelope(comp(), s, fb; envelope = :amplitude)[k,:]
    @test ea ≈ sqrt.(e) rtol = 1e-3
    @test !isapprox(ea, sqrt.(e); rtol = 1e-9)
    @test sum(ea)/length(ea) ≈ A rtol = 0.01

    @test_throws ErrorException power_envelope(comp(), s, fb; envelope = :hilbert)
end


@testset "decimation" begin
    fb = mbank()
    s = signal(am_tone(1.0, 0.5, 200.0, 7.0, 20.0), mfs)

    #max_mod_frq sets the decimation factor as floor(fs/(2.5*max_mod_frq)); fs_env overrides it
    for (mmf, M) in ((125, 10), (50, 25), (10, 125))
        tf = power_envelope(s, fb; max_mod_frq = mmf)
        @test tf.time[2] - tf.time[1] ≈ M/mfs
        @test power_envelope(s, fb; fs_env = mfs/M).components == tf.components
    end

    #An envelope rate that does not divide the signal rate is refused, with the neighbours named
    @test_throws ErrorException power_envelope(s, fb; fs_env = 300)
    @test_throws ErrorException power_envelope(s, fb; fs_env = 0)
    @test_throws ErrorException power_envelope(s, fb; max_mod_frq = 0)

    #The anti-alias filter has unit DC gain, so decimating does not move the envelope's mean
    und = power_envelope(comp(), s, fb; fs_env = mfs, trim = false)
    dec = power_envelope(comp(), s, fb; fs_env = mfs/25, trim = false)
    k = frqindex(200.0, fb.fcs)
    @test sum(dec[k,:])/size(dec,2) ≈ sum(und[k,:])/size(und,2) rtol = 1e-3

    #M = 1 is no decimation at all: no resampling filter is designed and no transient introduced
    @test size(und, 2) == length(s.x) - compensation_lead(gammatone_delay(fb; delay = :auto))
end


#The decimation is only legitimate because the resampler filters before it decimates. Modulate
#above the envelope Nyquist and no alias may appear at the frequency the fold would put it.
@testset "anti-aliasing" begin
    #One wide channel at 1 kHz: bw = 400 Hz passes the +/-100 Hz sidebands of the fast modulation
    fb = gammatone_filterbank(mfs, [1000.0]; bw = [400.0])
    fs_env = mfs/25                       #125 Hz, so the envelope Nyquist is 62.5 Hz

    fast = signal(am_tone(1.0, 0.5, 1000.0, 100.0, 20.0), mfs)   #100 Hz folds to |125-100| = 25
    slow = signal(am_tone(1.0, 0.5, 1000.0, 25.0, 20.0), mfs)    #the same rate, genuinely present

    pf = welch_psd(comp(), power_envelope(comp(), fast, fb; fs_env)[1,:], fs_env; frame_dur = 4.0)[1]
    ps, f = welch_psd(comp(), power_envelope(comp(), slow, fb; fs_env)[1,:], fs_env; frame_dur = 4.0)
    j = frqindex(25.0, f)

    @test ps[j] > 1e3*pf[j]                       #the genuine 25 Hz line dwarfs any alias
    @test pow2db(pf[j]/ps[j]) < -30
    @test argmax(ps) == j                         #and the genuine one is the largest peak
end


#A channel of bandwidth B can only show a modulation at f_m if both sidebands fc +/- f_m fit
#inside it, so the modulation only appears once B is around 2*f_m. No citation was found for this
#rule in the literature, so it is pinned here by measurement rather than asserted in a comment;
#it is what the bandwidth floor in SONAR's modspec rests on.
@testset "bandwidth sets the visible modulation rate" begin
    fc, fm, m = 500.0, 10.0, 0.5
    s = signal(am_tone(1.0, m, fc, fm, 20.0), mfs)

    depth = Float64[]
    for ratio in (0.5, 1.0, 2.0, 4.0)
        fb = gammatone_filterbank(mfs, [fc]; bw = [ratio*fm])
        e = power_envelope(comp(), s, fb; envelope = :amplitude, max_mod_frq = 50)[1,:]
        push!(depth, (maximum(e) - minimum(e))/(maximum(e) + minimum(e)))
    end

    @test issorted(depth)                 #wider channels always show more of the modulation
    @test depth[1] < 0.1*m                #B = f_m/2: essentially nothing gets through
    @test depth[4] > 0.7*m                #B = 4*f_m: most of the true depth is recovered
    @test depth[3] > 5*depth[1]           #the rule bites between B = f_m and B = 2*f_m
end


@testset "alignment" begin
    fb = mbank()
    s = signal(am_tone(1.0, 0.5, 200.0, 7.0, 20.0), mfs)
    d = gammatone_delay(fb; delay = :auto)
    lead = compensation_lead(d)

    #Compensation only shifts and phase-rotates, and the envelope discards phase, so slicing each
    #channel at its own offset must give exactly what compensating and trimming the complex
    #subbands gives - with no fabricated sample anywhere in it
    E = power_envelope(comp(), s, fb; fs_env = mfs, trim = false, align = d)
    Yc = compensate(gammatone_analysis(comp(), s, fb), d; trim = true)
    @test size(E) == size(Yc)
    @test E ≈ abs2.(Yc)
    @test all(E[:,1] .> 0)

    #align = :auto is the default and equals passing the object it builds
    @test power_envelope(comp(), s, fb; fs_env = mfs, trim = false) == E
    @test power_envelope(comp(), s, fb; fs_env = mfs, trim = false, align = :auto) == E

    #Without compensation nothing is dropped and the rows keep their own delays
    E0 = power_envelope(comp(), s, fb; fs_env = mfs, trim = false, align = nothing)
    @test size(E0, 2) == length(s.x)
    @test size(E0, 2) - size(E, 2) == lead

    #What alignment buys: a click peaks at the same instant in every channel with it, and at
    #visibly different instants without it
    click = zeros(4000); click[500] = 1.0
    sc = signal(click, mfs)
    pk_al = [argmax(view(power_envelope(comp(), sc, fb; fs_env = mfs, trim = false), k, :)) for k in eachindex(fb.fcs)]
    pk_no = [argmax(view(power_envelope(comp(), sc, fb; fs_env = mfs, trim = false, align = nothing), k, :)) for k in eachindex(fb.fcs)]
    @test length(unique(pk_al)) == 1
    @test length(unique(pk_no)) > 1
    @test maximum(pk_no) - minimum(pk_no) > 5

    #An explicit target too short to hold the slowest channels says so
    @test_logs (:warn, r"stay misaligned") power_envelope(comp(), s, fb; fs_env = mfs, align = 0.0005)
    @test_throws ErrorException power_envelope(comp(), s, fb; align = :nearest)
end


@testset "trimming and the time axis" begin
    fb = mbank()
    s = signal(am_tone(1.0, 0.5, 200.0, 7.0, 20.0), mfs)
    M = 10
    lead = compensation_lead(gammatone_delay(fb; delay = :auto))
    startup = startup_samples(fb)
    @test startup > lead                        #the start-up allowance is the binding one here

    tf = power_envelope(s, fb; max_mod_frq = 125)
    @test tf isa timefreq{Float64,Float64}
    @test tf.frqs == fb.fcs
    @test size(tf.components, 1) == length(fb.fcs)
    @test length(tf.time) == size(tf.components, 2)

    #The time axis keeps the signal's own origin: it starts where trimming left off, not at zero
    margin = TimeFrequencyAnalysis._resample_margin(Float64, M, 1, 60)
    @test tf.time[1] ≈ (startup + margin*M)/mfs
    @test tf.time[1] > 0
    @test tf.time[2] - tf.time[1] ≈ M/mfs

    #Untrimmed keeps everything the alignment left, and starts at the alignment lead
    raw = power_envelope(s, fb; max_mod_frq = 125, trim = false)
    @test raw.time[1] ≈ lead/mfs
    @test size(raw.components, 2) > size(tf.components, 2)
    @test occursin("Power Envelope", tf.title)
    @test occursin("Amplitude Envelope", power_envelope(s, fb; envelope = :amplitude).title)

    #A signal too short to survive trimming is an error that names the cause
    @test_throws ErrorException power_envelope(signal(randn(50), mfs), fb)
end


@testset "matrix form and precision" begin
    fb = mbank()
    s = signal(am_tone(1.0, 0.5, 200.0, 7.0, 20.0), mfs)
    d = gammatone_delay(fb; delay = :auto)

    #Given the complex subbands, the matrix form must reproduce the signal form exactly when the
    #latter is asked not to apply its filterbank-derived start-up trim
    Y = gammatone_analysis(comp(), s, fb)
    Em = power_envelope(comp(), Y, mfs; align = d, max_mod_frq = 125)
    Es = power_envelope(comp(), s, fb; align = d, max_mod_frq = 125, trim = false)
    #trim = false also removes the resample margin, so compare the matrix form without it too
    Em0 = power_envelope(comp(), Y, mfs; align = d, max_mod_frq = 125, trim = false)
    @test Em0 ≈ Es
    @test size(Em, 2) < size(Em0, 2)          #trim removes the resampler transients
    @test_throws ErrorException power_envelope(comp(), Y[1:3,:], mfs; align = d)

    #Precision follows the signal all the way through
    s32 = signal(Float32.(s.x), Float32(mfs))
    fb32 = gammatone_filterbank(Float32(mfs); fmin = 100, fmax = 400, base_frq = 200, filters_per_erb = 2, bw_scale = 2)
    tf32 = power_envelope(s32, fb32)
    @test tf32 isa timefreq{Float32,Float32}
    @test eltype(tf32.time) == Float32
    @test isapprox(tf32.components, Float32.(power_envelope(s, fb).components); rtol = 1e-3)
end
