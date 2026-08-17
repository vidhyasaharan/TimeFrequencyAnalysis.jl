#Visual verification companion to the modulation spectrum plan
#(dev docs/modulation-spectrum.md, "Build order"): one section per implementation step, each
#producing the plots that step's "Verify visually" item calls for. Run with Plots installed
#(Plots is not a package dependency): julia> include("scripts/modulation-verification.jl")
#Not part of the package or its test suite; by the final step this doubles as the source of
#figures for the docs page.
#
#Set TFA_VERIFY_FIGDIR to also write each figure to disk:
#    TFA_VERIFY_FIGDIR=/tmp/figs julia --project -e 'include("scripts/modulation-verification.jl")'

using TimeFrequencyAnalysis
using Plots
using Random
using Statistics

Random.seed!(2026)

const FIGDIR = get(ENV, "TFA_VERIFY_FIGDIR", "")
function keep(p, name)
    isempty(FIGDIR) || (mkpath(FIGDIR); savefig(p, joinpath(FIGDIR, name*".png")))
    return p
end

## ---------------------------------------------------------------------------------------
## Step 3 — welch_psd
## ---------------------------------------------------------------------------------------

#3a. Averaging is the whole point. A single periodogram of noise has a standard error of about
#100% of its own value HOWEVER LONG the record: lengthening the frame buys frequency resolution,
#not reliability. Only averaging F frames reduces the scatter, by 1/sqrt(F). Each curve below
#uses the same 200 s of white noise, cut into F frames at 50% overlap, so they differ only in
#how many estimates were averaged - and the mean must stay on the true level throughout.
step3a = let
    fs = 1000.0
    N = 200_000
    σ = 2.0
    w = σ.*randn(N)
    truth = 2σ^2/fs                       #two-sided sigma^2/fs folded onto positive frequencies

    p = plot(; xlabel = "Frequency (Hz)", ylabel = "PSD (dB re 1/Hz)",
             title = "welch_psd: averaging F frames cuts the scatter by 1/√F",
             legend = :bottomright, size = (900, 520), ylims = (-60, -5))
    for (F, col) in zip((1, 4, 16, 64), (:grey70, :orange, :seagreen, :navy))
        dur = 2*(N/fs)/(F+1)              #50% overlap: F frames need frame_dur = 2*total/(F+1)
        sp = welch_psd(w, fs; frame_dur = dur, scaling = :density)
        b = sp.components[5:end-5]
        m = sum(b)/length(b)
        sd = sqrt(sum(abs2, b .- m)/length(b))/m
        plot!(p, sp.frqs, pow2db.(sp.components);
              label = "F = $F  (scatter $(round(100sd, digits = 0))%)", lc = col, lw = F == 1 ? 0.6 : 1.2)
    end
    hline!(p, [pow2db(truth)]; lc = :red, ls = :dash, lw = 2, label = "true 2σ²/fs")
    keep(p, "step3a-variance-reduction")
end
display(step3a)

#3b. The two scaling conventions answer different questions, and the difference only shows up
#when the frame length changes. :spectrum holds a LINE fixed at A^2/2 while the noise floor
#drops as the frames lengthen (each bin collects less noise); :density holds the NOISE FLOOR
#fixed while the line climbs (a line's power concentrates into an ever narrower bin). Reading
#a line off a density, or a noise floor off a spectrum, is the classic PSD mistake - here it is
#visible rather than merely asserted.
step3b = let
    fs = 1000.0
    N = 200_000
    t = (0:N-1)./fs
    A = 1.7
    x = A.*cos.(2π*50 .*t) .+ 0.3.*randn(N)

    ps = plot(; ylabel = "Power (dB re 1)", title = "scaling = :spectrum — the LINE is invariant",
              legend = :topright, xlims = (30, 70))
    pd = plot(; xlabel = "Frequency (Hz)", ylabel = "PSD (dB re 1/Hz)",
              title = "scaling = :density — the NOISE FLOOR is invariant",
              legend = :topright, xlims = (30, 70))
    for (fd, col) in zip((0.5, 1.0, 2.0), (:orange, :seagreen, :navy))
        a = welch_psd(x, fs; frame_dur = fd, scaling = :spectrum)
        b = welch_psd(x, fs; frame_dur = fd, scaling = :density)
        plot!(ps, a.frqs, pow2db.(a.components); label = "frame $fd s", lc = col)
        plot!(pd, b.frqs, pow2db.(b.components); label = "frame $fd s", lc = col)
    end
    hline!(ps, [pow2db(A^2/2)]; lc = :red, ls = :dash, lw = 2, label = "A²/2")
    hline!(pd, [pow2db(2*0.3^2/fs)]; lc = :red, ls = :dash, lw = 2, label = "2σ²/fs")
    keep(plot(ps, pd; layout = (2,1), size = (900, 640)), "step3b-scaling-conventions")
end
display(step3b)

#3c. The four dc modes, on a strictly positive envelope of mean 2.5 modulated 40% at 3 Hz.
#Two independent operations, and they do NOT behave alike:
#  * dividing by the mean rescales EVERY bin by exactly 1/mu^2, whatever the window, so the
#    :divide/:keep and :index/:remove ratios are flat lines at 1/mu^2 (lower panel);
#  * subtracting the mean reaches only the DC bin under a rectangular window, but a tapered
#    window smears the frame's constant over its mainlobe and sidelobes. That is why :keep sits
#    tens of dB above :remove in the bins either side of DC - exactly where shaft lines live.
#Removing the mean BEFORE windowing is what suppresses that leakage.
step3c = let
    fs = 100.0
    N = 6000
    t = (0:N-1)./fs
    μ = 2.5
    m = 0.4
    env = μ.*(1 .+ m.*cos.(2π*3 .*t))
    kw = (frame_dur = 4.0, scaling = :spectrum)

    P = Dict(d => welch_psd(comp(), env, fs; dc = d, kw...)[1] for d in (:keep,:remove,:divide,:index))
    f = welch_psd(comp(), env, fs; kw...)[2]

    ptop = plot(; ylabel = "Power (dB)", title = "welch_psd dc modes (Hann window, μ = $μ, 40% at 3 Hz)",
                legend = :topright, xlims = (0, 12), ylims = (-130, 25))
    for (d, col) in zip((:keep,:remove,:divide,:index), (:black,:crimson,:steelblue,:seagreen))
        plot!(ptop, f, pow2db.(P[d]); label = ":$d", lc = col, lw = 1.6)
    end
    vline!(ptop, [3.0]; lc = :grey, ls = :dot, label = "modulation rate")
    annotate!(ptop, 4.2, -15,
              text("black above red below ~2.5 Hz = the mean's\nleakage that :remove suppresses", 7, :left, :black))

    #Division is exact: both ratios are flat at 1/mu^2 across the whole axis
    pmid = plot(; ylabel = "ratio", legend = :topright, xlims = (0, 12), ylims = (0, 0.32))
    plot!(pmid, f, P[:divide]./max.(P[:keep], eps());   label = ":divide / :keep",   lc = :steelblue, lw = 1.6)
    plot!(pmid, f, P[:index]./max.(P[:remove], eps());  label = ":index / :remove",  lc = :seagreen, lw = 1.6, ls = :dash)
    hline!(pmid, [1/μ^2]; lc = :red, ls = :dash, lw = 2, label = "1/μ² = $(round(1/μ^2, digits = 3))")

    #Subtraction is window dependent: rect confines it to DC, Hann does not
    pbot = plot(; xlabel = "Frequency (Hz)", ylabel = "Power (dB)",
                title = "subtraction is DC-only for rect, not for a taper", legend = :topright,
                xlims = (0, 12), ylims = (-220, 20))
    for (wt, ls) in (("rect", :solid), ("hanning", :dash))
        for (d, col) in ((:keep, :black), (:remove, :crimson))
            Q = welch_psd(comp(), env, fs; dc = d, wtype = wt, kw...)[1]
            plot!(pbot, f, pow2db.(Q); label = ":$d, $wt", lc = col, ls = ls, lw = 1.4)
        end
    end
    keep(plot(ptop, pmid, pbot; layout = (3,1), size = (900, 900)), "step3c-dc-modes")
end
display(step3c)

#3d. Why 3c looks the way it does, and the framing welch_psd actually used. Left: the windows
#and their transforms - the Hann mainlobe is four bins wide against the rectangular window's
#two, which is the leakage seen above, bought in exchange for sidelobes 30 dB lower. Right: the
#requested framing against the framing that ran. welch_psd takes its geometry from
#framed_signal/view_frame and averages only the num_signal_frames whole frames, so the shift is
#exactly the one asked for and no zero padded tail frame dilutes the average.
step3d = let
    L = 256
    pw = plot(; xlabel = "Sample", ylabel = "w[n]", title = "Windows", legend = :topright)
    pr = plot(; xlabel = "Bins from centre", ylabel = "dB", title = "Window transforms",
              legend = :topright, xlims = (-8, 8), ylims = (-90, 5))
    for (wt, col) in (("rect", :black), ("hanning", :seagreen), ("hamming", :orange))
        w = window(L; wtype = wt)
        enbw = L*sum(abs2, w)/sum(w)^2
        plot!(pw, 0:L-1, w; label = "$wt (ENBW $(round(enbw, digits = 2)) bins)", lc = col)
        W = abs.(dft(w; ndft = 8L, wtype = "rect"))
        plot!(pr, (0:length(W)-1)./8, amp2db.(W./maximum(W)); label = wt, lc = col)
        plot!(pr, -(0:length(W)-1)./8, amp2db.(W./maximum(W)); label = "", lc = col)
    end

    fs = 100.0
    s = signal(randn(6050), fs)
    rows = String[]
    for (fd, fsh) in ((4.0, 2.0), (4.0, 1.0), (3.3, 1.65))
        fr = framed_signal(s, fd, fsh)
        starts = [(i-1)*fr.frame_shift + 1 for i in 1:fr.num_signal_frames]
        push!(rows, "$(fd)s / $(fsh)s → shift $(fr.frame_shift) samples " *
                    "(asked $(round(Int, fsh*fs))), $(fr.num_signal_frames) whole of $(fr.num_frames), " *
                    "unique gaps $(unique(diff(starts)))")
    end
    pt = plot(; framestyle = :none, title = "Framing requested vs framing used", legend = false)
    for (i, r) in enumerate(rows)
        annotate!(pt, 0.0, 1.0 - 0.18i, text(r, 7, :left, :black))
    end
    xlims!(pt, 0, 1); ylims!(pt, 0, 1)

    keep(plot(pw, pr, pt; layout = @layout([a b; c{0.25h}]), size = (1000, 620)), "step3d-window-and-framing")
end
display(step3d)

#3e. What averaging does and does not buy for DETECTION. The received wisdom is "average to pull
#a tone out of the noise", and for a STABLE tone in WHITE noise from a FIXED record that is
#simply wrong: one long transform is the matched filter and wins outright. With
#scaling = :spectrum the line stays at A^2/2 while the mean floor falls as 1/L, so halving the
#number of frames doubles the frame length and buys 3 dB of processing gain - more than the
#tighter detection threshold that averaging allows can repay. Going F = 1 -> F = 64 below gives
#up 15.1 dB of processing gain and recovers only 7.0 dB from the threshold, a net 8.1 dB loss;
#the panels measure 13.7 -> 5.5 dB, i.e. 8.2 dB. The caption on an earlier draft of this figure
#claimed the opposite, which is exactly what a visual check is for.
#
#What averaging DOES buy is the floor itself. At F = 1 the floor scatters by about 100% of its
#own value, so it is useless as a level estimate, as a display scale, or as a threshold that
#means anything; by F = 64 it is a measurement good to ~13%. That, plus robustness when the line
#wanders (3f) or the scene is not stationary, is why modulation_spectrum averages frames instead
#of transforming the whole span at once.
step3e_detect = let
    fs = 1000.0
    N = 100_000                                   #100 s of data, fixed: this is what is on offer
    t = (0:N-1)./fs
    f0 = 50.0
    σ = 1.0
    A = 0.1                                       #per-sample SNR of A^2/2/σ^2 = -23 dB
    x = A.*cos.(2π*f0.*t) .+ σ.*randn(N)

    p = plot(; layout = (3,1), size = (900, 900), legend = :topright)
    for (row, F) in enumerate((1, 8, 64))
        dur = 2*(N/fs)/(F+1)
        sp = welch_psd(x, fs; frame_dur = dur, scaling = :spectrum)
        band = abs.(sp.frqs .- f0) .> 2.0         #noise-only bins, tone excluded
        noise = sp.components[band]
        m = mean(noise)
        scatter_pct = 100*std(noise)/m            #how usable the floor is as an estimate
        thr = quantile(noise, 0.999)              #empirical 1-in-1000 false alarm threshold
        pk, ipk = findmax(sp.components)
        margin = pow2db(pk/thr)
        plot!(p, sp.frqs, pow2db.(sp.components); subplot = row, lc = :steelblue, lw = 0.7,
              label = "F = $F", xlims = (0, 500), ylabel = "Power (dB)",
              title = "F = $F, frame $(round(dur, digits = 1)) s — detection margin " *
                      "$(round(margin, digits = 1)) dB, floor scatter $(round(Int, scatter_pct))%",
              titlefontsize = 9)
        hline!(p, [pow2db(thr)]; subplot = row, lc = :red, ls = :dash, lw = 1.5,
               label = "1-in-1000 false alarm")
        hline!(p, [pow2db(m)]; subplot = row, lc = :grey, ls = :dot, lw = 1.2, label = "mean floor")
        scatter!(p, [sp.frqs[ipk]], [pow2db(pk)]; subplot = row, mc = :orange, ms = 5,
                 markerstrokewidth = 0, label = "peak")
    end
    xlabel!(p[3], "Frequency (Hz)")
    annotate!(p[1], 250, -32,
              text("margin FALLS with averaging here; what improves is the floor", 8, :center, :black))
    keep(p, "step3e-detection-in-noise")
end
display(step3e_detect)

#3f. How long should the frames be? For a PERFECTLY stable tone, longer is monotonically better:
#with scaling = :spectrum the line stays at A^2/2 while the noise floor falls as 1/L, so the
#processing gain grows with frame length faster than the loss of averages hurts. Real machinery
#lines are not stable - a shaft rate wanders - and once the resolution is finer than the line's
#own wander the tone smears across bins and lengthening the frame only costs averages. That puts
#a plateau in the curve, and it is the reason modulation_spectrum defaults to a 10 s
#frame rather than simply the longest one that fits.
step3f_tradeoff = let
    fs = 1000.0
    N = 200_000
    t = (0:N-1)./fs
    f0 = 50.0
    σ = 1.0
    A = 0.1

    stable = A.*cos.(2π*f0.*t)
    #A line that wanders +/- 0.4 Hz over a 25 s cycle: integrate the instantaneous frequency
    finst = f0 .+ 0.4.*sin.(2π.*t./25)
    drifting = A.*cos.(2π.*cumsum(finst)./fs)

    durs = exp10.(range(log10(0.25), log10(40.0); length = 22))
    margins = Dict(k => Float64[] for k in ("stable","wandering"))
    for (name, tone) in (("stable", stable), ("wandering", drifting))
        y = tone .+ σ.*randn(N)
        for d in durs
            sp = welch_psd(y, fs; frame_dur = d, scaling = :spectrum)
            band = abs.(sp.frqs .- f0) .> 3.0
            sorted = sort(sp.components[band])
            thr = sorted[max(1, round(Int, 0.999*length(sorted)))]
            near = abs.(sp.frqs .- f0) .<= 3.0
            push!(margins[name], pow2db(maximum(sp.components[near])/thr))
        end
    end

    p = plot(durs, margins["stable"]; xscale = :log10, lc = :navy, lw = 2, marker = :circle, ms = 3,
             label = "stable line", xlabel = "Frame duration (s)",
             ylabel = "Detection margin (dB above 1-in-1000 threshold)",
             title = "Frame length: longer always wins for a stable line, not for a wandering one",
             titlefontsize = 10, legend = :bottomright, size = (900, 520))
    plot!(p, durs, margins["wandering"]; lc = :crimson, lw = 2, marker = :circle, ms = 3,
          label = "line wandering ±0.4 Hz")
    ibest = argmax(margins["wandering"])
    vline!(p, [durs[ibest]]; lc = :crimson, ls = :dash,
           label = "best for wandering ≈ $(round(durs[ibest], digits = 1)) s")
    hline!(p, [0.0]; lc = :grey, ls = :dot, label = "detection threshold")
    keep(p, "step3f-frame-length-tradeoff")
end
display(step3f_tradeoff)

#3g. What the averaging costs. Two tones 1.5 Hz apart in noise: a long frame resolves them but
#leaves a floor too ragged to threshold, a short frame gives a floor you can threshold but merges
#them into one line. welch_psd does not escape the resolution/variance trade - it just puts the
#dial where the caller can reach it, via frame_dur.
step3g_resolution = let
    fs = 1000.0
    N = 200_000
    t = (0:N-1)./fs
    x = 0.12.*cos.(2π*50.0.*t) .+ 0.12.*cos.(2π*51.5.*t) .+ randn(N)

    p = plot(; layout = (2,1), size = (900, 620), legend = :topright)
    for (row, (d, note)) in enumerate(((20.0, "resolved, but the floor is ragged"),
                                       (1.0,  "smooth floor, but the two lines have merged")))
        sp = welch_psd(x, fs; frame_dur = d, scaling = :spectrum)
        F = framed_signal(signal(x, fs), d, d/2).num_signal_frames
        plot!(p, sp.frqs, pow2db.(sp.components); subplot = row, lc = :steelblue, lw = 1.2,
              xlims = (45, 57), label = "frame $d s  (F = $F, Δf = $(round(1/d, digits = 2)) Hz)",
              title = "$note", titlefontsize = 9, ylabel = "Power (dB)")
        vline!(p, [50.0, 51.5]; subplot = row, lc = :red, ls = :dot, label = "true lines")
    end
    xlabel!(p[2], "Frequency (Hz)")
    keep(p, "step3g-resolution-variance-tradeoff")
end
display(step3g_resolution)
