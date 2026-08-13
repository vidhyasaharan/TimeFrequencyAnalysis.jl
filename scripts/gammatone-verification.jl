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
    t, h = impulse_response(gf)
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
    t, h = impulse_response(gf; dur = N/fs)
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
