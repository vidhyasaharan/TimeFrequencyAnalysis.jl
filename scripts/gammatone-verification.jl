#Visual verification companion to the gammatone filterbank plan
#(dev docs/gammatone-filterbank.md §4): one section per implementation step, each producing
#the plots that step's "Verify visually" item calls for. Run with Plots installed (Plots is
#not a package dependency): julia> include("scripts/gammatone-verification.jl")
#Not part of the package or its test suite; by the final step this doubles as the source of
#figures for the docs page.

using TimeFrequencyAnalysis
using Plots

## ---------------------------------------------------------------------------------------
## Step 1 — generalised filter_coefs: a complex one-pole filter has a one-sided response
## ---------------------------------------------------------------------------------------
#A single complex pole ã = λ e^{iβ} resonates only at +fp: the magnitude response has one
#asymmetric peak of gain 1/(1-λ) at the pole frequency and no mirror peak at fs-fp (i.e.
#-fp), which no real-coefficient filter can do. The phase now follows the standard z⁻¹
#convention (e.g. a pure delay has phase -2πf/fs), so phase plots and group delays read
#correctly for the gammatone work that builds on this.
step1 = let
    fs = 8000.0
    λ = 0.9
    fp = 1000.0
    F = filter_coefs([1.0], [1.0, -λ*exp(im*2π*fp/fs)])
    f = collect(0.0:5.0:fs-5.0)          #full [0, fs) to show the absence of a mirror peak
    resp = filter_resp(F, f, fs)

    pmag = plot(f, amp2db.(abs.(resp));
                ylabel = "Magnitude (dB)", label = "",
                title = "Complex one-pole, fp = 1 kHz, λ = 0.9: one-sided resonance")
    vline!(pmag, [fp, fs - fp]; ls = :dash, label = "fp and fs-fp (mirror)")

    pph = plot(f, angle.(resp);
               xlabel = "Frequency (Hz)", ylabel = "Phase (rad)", label = "")

    plot(pmag, pph; layout = (2,1), size = (800, 550))
end
display(step1)

## ---------------------------------------------------------------------------------------
## Step 2 — gammatone_filter design: impulse response, magnitude and phase responses
## ---------------------------------------------------------------------------------------
#The 4th-order gammatone at 1 kHz (fs = 16 kHz): the impulse response is a tone burst under
#a gamma-shaped envelope peaking at (γ-1)/(2πb) ≈ 3.5 ms; the magnitude response peaks at
#fc with gain exactly 2 (+6.02 dB); the phase falls through the passband.
step2a = let
    fs = 16000.0
    gf = gammatone_filter(fs, 1000.0)
    t, h = gammatone_impulse_response(gf)
    p1 = plot(1000 .* t, [real.(h) abs.(h)];
              label = ["real part" "envelope"], xlabel = "Time (ms)",
              title = "Impulse response, fc = 1 kHz, order 4, bw = ERB(fc) ≈ 132.6 Hz")
    f = collect(0.0:2.0:8000.0)
    resp = filter_resp(gf, f)
    p2 = plot(f, amp2db.(abs.(resp)); label = "", ylabel = "Magnitude (dB)", ylims = (-80, 10))
    vline!(p2, [1000.0]; ls = :dash, label = "fc (peak +6.02 dB)")
    p3 = plot(f, angle.(resp); label = "", xlabel = "Frequency (Hz)", ylabel = "Phase (rad)")
    plot(p1, p2, p3; layout = (3,1), size = (820, 850))
end
display(step2a)

#Design variants: Route A (ERB-matching, the default) vs Route B (edges at -3 dB) for the
#same nominal bandwidth; the bw_scale family; explicit bw overriding the ERB default.
step2b = let
    fs = 16000.0
    f = collect(200.0:2.0:2600.0)
    B = erb(1000.0)

    gA = gammatone_filter(fs, 1000.0)
    gB = gammatone_filter(fs, 1000.0, B, 3.0)
    p1 = plot(f, amp2db.(filter_magresp(gA, f)); label = "Route A: ERB(fc) as filter ERB",
              ylabel = "Magnitude (dB)", ylims = (-60, 10),
              title = "Same bandwidth value, two conventions")
    plot!(p1, f, amp2db.(filter_magresp(gB, f)); label = "Route B: ERB(fc) as -3 dB width")
    scatter!(p1, [1000 - B/2, 1000 + B/2], fill(amp2db(2) - 3, 2);
             label = "fc ± bw/2 at -3 dB (Route B)", marker = :circle)

    p2 = plot(; ylabel = "Magnitude (dB)", ylims = (-60, 10), title = "bw_scale ∈ {0.5, 1, 2}")
    for s in (0.5, 1.0, 2.0)
        g = gammatone_filter(fs, 1000.0; bw_scale = s)
        plot!(p2, f, amp2db.(filter_magresp(g, f)); label = "bw_scale = $s")
    end

    p3 = plot(; xlabel = "Frequency (Hz)", ylabel = "Magnitude (dB)", ylims = (-60, 10),
              title = "Explicit bw ∈ {80, ERB(fc) ≈ 132.6, 300} Hz")
    for bwv in (80.0, round(erb(1000.0); digits = 1), 300.0)
        g = gammatone_filter(fs, 1000.0; bw = bwv)
        plot!(p3, f, amp2db.(filter_magresp(g, f)); label = "bw = $bwv Hz")
    end

    plot(p1, p2, p3; layout = (3,1), size = (820, 850))
end
display(step2b)

## ---------------------------------------------------------------------------------------
## Step 3 — filtering kernel: recursion vs closed form, and steady-state behaviour
## ---------------------------------------------------------------------------------------
#Filtering a unit impulse with the recursive kernel must land on the closed-form impulse
#response exactly (they are independent computations of the same thing); a tone at fc must
#settle to envelope 1 (gain 2 halved across ±fc); and on a two-tone-in-noise signal the
#channel envelopes must sit at their components' amplitudes.
step3 = let
    fs = 16000.0
    gf = gammatone_filter(fs, 1000.0)

    #Kernel impulse output vs closed form (markers subsampled so the overlay stays legible)
    N = 400
    e = zeros(N); e[1] = 1.0
    y = gammatone_filt(gf, e)
    t, h = gammatone_impulse_response(gf; dur = N/fs)
    maxerr = maximum(abs.(y .- h))
    p1 = plot(1000 .* t, real.(h); label = "closed form (real part)",
              title = "Kernel vs closed form: max |difference| = $(round(maxerr, sigdigits=3))")
    scatter!(p1, 1000 .* t[1:8:end], real.(y[1:8:end]); label = "kernel output",
             marker = :circle, ms = 2.5, msw = 0)
    plot!(p1, 1000 .* t, abs.(h); label = "envelope", ls = :dash)

    #Tone at fc: envelope settles to 1 after the ~group-delay transient
    n = 0:Int(round(0.05*fs))-1
    x = cos.(2π*1000.0.*n./fs)
    env = abs.(gammatone_filt(gf, x))
    p2 = plot(1000 .* n./fs, env; label = "envelope of filtered 1 kHz tone",
              xlabel = "Time (ms)")
    hline!(p2, [1.0]; ls = :dash, label = "unit amplitude")

    #Two tones in noise (the test-suite fixture recipe, fs = 8 kHz): each channel's
    #envelope hovers at its component's amplitude
    fs8 = 8000.0
    n8 = 0:Int(round(0.5*fs8))-1
    x8 = cos.(2π*440.0.*n8./fs8) .+ 0.5.*cos.(2π*1000.0.*n8./fs8) .+ 0.3.*randn(length(n8))
    e440 = abs.(gammatone_filt(gammatone_filter(fs8, 440.0), x8))
    e1k = abs.(gammatone_filt(gammatone_filter(fs8, 1000.0), x8))
    p3 = plot(1000 .* n8./fs8, [e440 e1k];
              label = ["channel at 440 Hz" "channel at 1 kHz"],
              xlabel = "Time (ms)", title = "Two tones (1.0 at 440 Hz, 0.5 at 1 kHz) in noise")
    hline!(p3, [1.0, 0.5]; ls = :dash, label = "")

    plot(p1, p2, p3; layout = (3,1), size = (820, 850))
end
display(step3)

## ---------------------------------------------------------------------------------------
## Step 4 — filterbank and analyses: bank responses and the first cochleagrams
## ---------------------------------------------------------------------------------------
#The default 30-channel bank (70 Hz – 6.7 kHz, anchored at 1 kHz): magnitude responses
#widen in proportion to the ERB as centre frequency rises, while a constant-bw bank keeps
#every channel the same width. The cochleagram of a linear chirp then shows per-sample
#time resolution on an ERB frequency axis, next to the framewise, uniform-resolution
#specgram of the same signal.
step4 = let
    fs = 16000.0
    fb = gammatone_filterbank(fs)
    f = collect(0.0:4.0:8000.0)

    p1 = plot(; ylims = (-50, 10), ylabel = "Magnitude (dB)", legend = false,
              title = "Default bank: 30 channels, bw = ERB(fc)")
    for gf in fb.filters
        plot!(p1, f, amp2db.(filter_magresp(gf, f)))
    end

    fbc = gammatone_filterbank(fs; bw = 150.0)
    p2 = plot(; ylims = (-50, 10), xlabel = "Frequency (Hz)", ylabel = "Magnitude (dB)",
              legend = false, title = "Constant-bandwidth bank: bw = 150 Hz everywhere")
    for gf in fbc.filters
        plot!(p2, f, amp2db.(filter_magresp(gf, f)))
    end

    #Linear chirp, 100 Hz → 7 kHz over one second
    t = (0:Int(fs)-1)./fs
    chirp = sin.(2π.*(100.0.*t .+ 3450.0.*t.^2))
    s = signal(chirp, fs)
    #clims clip the eps-floor of quiet channels (-300 dB) so the useful range is visible
    p3 = plot(amp2db(gammatone_cochleagram(s)); clims = (-80, 0))
    p4 = plot(amp2db(specgram(framed_signal(s))); clims = (-45, 35))

    plot(p1, p2, p3, p4; layout = @layout([a; b; c d]), size = (950, 1000))
end
display(step4)

## ---------------------------------------------------------------------------------------
## Step 5 — delay helpers: group delay, envelope-peak delay, and the unaligned sum
## ---------------------------------------------------------------------------------------
#Gammatone channels are slow at the bottom of the bank and fast at the top (both delays
#scale as 1/b). The group delay at fc and the envelope peak are different instants — the
#envelope peaks earlier — and the per-channel spread is exactly why the raw summed bank
#response is a mess until channels are delay- and phase-aligned (step 6).
step5 = let
    fs = 16000.0
    fb = gammatone_filterbank(fs)

    #Group-delay curves for a low, mid and high channel (dots: closed form γλ/(1-λ) at fc)
    p1 = plot(; xlabel = "Frequency (Hz)", ylabel = "Group delay (ms)",
              title = "Group delay of a low, mid and high channel")
    for fc0 in (250.0, 1000.0, 4000.0)
        gf = fb.filters[frqindex(fc0, fb.fcs)]
        fgrid = collect(range(max(50.0, gf.fc - 2*gf.bw), gf.fc + 2*gf.bw; length = 400))
        plot!(p1, fgrid, 1000 .* group_delay(gf, fgrid)./fs; label = "fc ≈ $(round(Int, gf.fc)) Hz")
        scatter!(p1, [gf.fc], [1000*group_delay(gf)/fs]; label = "", marker = :circle, ms = 4)
    end

    #Envelope-peak delay across the bank, against the analogue (order-1)/(2πb) curve
    tpk = [1000*envelope_delay(g)/fs for g in fb.filters]
    tpk_analog = [1000*(g.order - 1)/(2π*g.b) for g in fb.filters]
    p2 = plot(fb.fcs, tpk_analog; xscale = :log10, label = "(order-1)/(2πb) (analogue)",
              xlabel = "Centre frequency (Hz)", ylabel = "Envelope peak (ms)",
              title = "Envelope-peak delay across the bank")
    scatter!(p2, fb.fcs, tpk; label = "measured argmax of |h|", ms = 3, msw = 0)

    #The uncompensated summed response: deep interference ripple — what the delay/phase
    #alignment (step 6) and the mixer (v2) exist to fix
    fgrid = collect(50.0:2.0:7950.0)
    p3 = plot(fgrid, amp2db.(abs.(summed_resp(fb, fgrid)));
              label = "", xlabel = "Frequency (Hz)", ylabel = "Magnitude (dB)",
              title = "Summed bank response without alignment")

    plot(p1, p2, p3; layout = (3,1), size = (850, 900))
end
display(step5)

## ---------------------------------------------------------------------------------------
## Step 6 — delay and phase compensation: the classic before/after alignment figure
## ---------------------------------------------------------------------------------------
#Left: raw per-channel impulse responses (real parts, each normalised) — the envelope
#peaks drift later and later towards the bottom of the bank. Right: after gammatone_delay
#compensation the envelopes and the fine structure line up at the 4 ms target (the slowest
#channels cannot be advanced and stay put). Below: the summed real parts become a pulse at
#4 ms. In the frequency domain what alignment buys is PHASE coherence — the compensated
#sum is a near-pure nd-sample delay where channels are aligned — while the magnitude keeps
#a picket-fence ripple (the incoherent raw sum is actually smoother in magnitude); that
#residue is the v2 mixer's job.
step6 = let
    fs = 16000.0
    fb = gammatone_filterbank(fs)
    d = gammatone_delay(fb; delay = 0.004)
    N = Int(round(0.03*fs))
    e = zeros(N); e[1] = 1.0
    Y = gammatone_filt(fb, e)
    Yc = compensate(Y, d)
    t = 1000 .* (0:N-1)./fs

    function stackplot(M, ttl)
        p = plot(; title = ttl, xlabel = "Time (ms)", ylabel = "Channel",
                 legend = false)
        for k in axes(M, 1)
            r = view(M, k, :)
            plot!(p, t, k .+ 0.6.*real.(r)./maximum(abs.(r)); c = :steelblue, lw = 0.8)
        end
        vline!(p, [1000*d.nd/fs]; ls = :dash, c = :black)
        p
    end
    p1 = stackplot(Y, "Impulse responses, raw")
    p2 = stackplot(Yc, "After delay + phase compensation")

    p3 = plot(t, [vec(sum(real.(Y); dims = 1)) vec(sum(real.(Yc); dims = 1))];
              label = ["raw sum" "compensated sum"], xlabel = "Time (ms)",
              title = "Σ of channel real parts: an impulse re-emerges at 4 ms")
    vline!(p3, [1000*d.nd/fs]; ls = :dash, c = :black, label = "")

    fgrid = collect(60.0:2.0:7900.0)
    θg = 2π.*fgrid./fs
    sr_raw = summed_resp(fb, fgrid)
    sr_cmp = summed_resp(fb, fgrid; delay = d)
    p4 = plot(fgrid, amp2db.(abs.(sr_raw));
              label = "raw", xlabel = "Frequency (Hz)", ylabel = "Magnitude (dB)",
              title = "Summed response magnitude")
    plot!(p4, fgrid, amp2db.(abs.(sr_cmp)); label = "delay compensated")
    p5 = plot(fgrid, angle.(sr_raw.*cis.(θg.*d.nd));
              label = "raw", xlabel = "Frequency (Hz)", ylabel = "Phase (rad)",
              title = "Summed-response phase after removing the 4 ms delay")
    plot!(p5, fgrid, angle.(sr_cmp.*cis.(θg.*d.nd)); label = "delay compensated", lw = 2)

    plot(p1, p2, p3, p4, p5; layout = @layout([a b; c; d e]), size = (1050, 1150))
end
display(step6)
