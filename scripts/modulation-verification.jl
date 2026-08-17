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


## ---------------------------------------------------------------------------------------
## Step 4 — modfreq
## ---------------------------------------------------------------------------------------

#4a. Orientation. A modfreq stores one ROW per carrier and one COLUMN per modulation rate, so
#the image must have carrier on y and modulation rate on x. A transposed heatmap still looks
#like a plausible modulation spectrum, which is exactly why it has to be checked against a
#pattern that is not symmetric: here component[k,j] = 10k + j, so the value rises slowly down
#the carrier axis and quickly along the modulation axis, and the corner values are all distinct.
#(The plot recipe lands at step 7; until then this draws the component matrix directly, which is
#what the recipe will wrap.)
step4a = let
    ncar, nmod = 6, 9
    C = [10.0k + j for k in 1:ncar, j in 1:nmod]
    fcs = [80.0, 120.0, 180.0, 270.0, 420.0, 640.0]          #six carriers, one per row
    mfr = collect(range(0.0, 40.0; length = nmod))
    mf = modfreq(C, fcs, mfr, 2 .*erb.(fcs), 312.5, 10.0, 11, "orientation check: 10k + j")

    p = heatmap(mf.components; size = (900, 460), c = :viridis,
                xticks = (1:nmod, string.(round.(Int, mf.mod_frqs))),
                yticks = (1:ncar, string.(round.(Int, mf.fcs))),
                xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                title = mf.title)
    for k in 1:ncar, j in 1:nmod
        annotate!(p, j, k, text(string(round(Int, C[k,j])), 6, :white))
    end
    #The corners must read 11 at (lowest carrier, lowest rate) and 69 at (highest, highest)
    keep(p, "step4a-orientation")
end
display(step4a)

#4b. Ticks on a log carrier axis. The carrier grid is spaced on the ERB rate scale, so it is far
#from uniform in Hz: drawn in index space with generate_ticks the rows stay evenly spaced on
#screen while the labels report the true frequencies, which is how the timefreq recipe already
#handles a nonuniform frequency axis. Plotting against the Hz values directly instead would
#crush the 16 carriers of band X into the bottom eighth of the image.
step4b = let
    fcs = erbfreq_array(; fmin = 77.0, fmax = 1232.0, base_frq = 308.0, filters_per_erb = 8)
    ncar = length(fcs)
    mfr = collect(range(0.0, 60.0; length = 120))
    #A ridge whose modulation rate rises linearly with CARRIER INDEX. In index space that is a
    #straight diagonal; against Hz the ERB grid bunches the low carriers together, so the lower
    #half of the diagonal collapses into the bottom of the image.
    C = [exp(-((mfr[j] - (3.0 + 45.0*(k-1)/(ncar-1)))/0.8)^2) + 0.01 for k in 1:ncar, j in eachindex(mfr)]
    mf = modfreq(C, fcs, mfr, 2 .*erb.(fcs), 312.5, 10.0, 11, "")

    ytk = TimeFrequencyAnalysis.generate_ticks(mf.fcs, 8)
    xtk = TimeFrequencyAnalysis.generate_ticks(mf.mod_frqs, 6)
    pidx = heatmap(mf.components; c = :viridis, xticks = xtk, yticks = ytk, colorbar = false,
                   xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                   title = "index space: $(ncar) evenly spaced rows", titlefontsize = 10,
                   left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    phz = heatmap(mf.mod_frqs, mf.fcs, mf.components; c = :viridis, colorbar = false,
                  xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                  title = "against Hz: the lower 60 rows collapse", titlefontsize = 10,
                  left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    keep(plot(pidx, phz; layout = (1,2), size = (1100, 480)), "step4b-carrier-axis-ticks")
end
display(step4b)

#4c. The element-wise operations must return a modfreq with every axis and every piece of
#provenance intact, not a bare matrix: a dB view of a modulation spectrum is still a modulation
#spectrum and must keep plotting like one. pow2db is 10log10 and amp2db 20log10, so on the same
#data one is exactly half the other - the right one to use here is pow2db, because
#modulation_spectrum returns power.
step4c = let
    fcs = erbfreq_array(; fmin = 77.0, fmax = 1232.0, base_frq = 308.0, filters_per_erb = 4)
    ncar = length(fcs)
    mfr = collect(range(0.0, 60.0; length = 50))
    C = [exp(-((mfr[j]-7.0)/1.0)^2)*(0.2 + k/ncar) + 0.001 for k in 1:ncar, j in 1:length(mfr)]
    mf = modfreq(C, fcs, mfr, 2 .*erb.(fcs), 312.5, 10.0, 11, "linear power")
    db = pow2db(mf)

    ytk = TimeFrequencyAnalysis.generate_ticks(mf.fcs, 6)
    xtk = TimeFrequencyAnalysis.generate_ticks(mf.mod_frqs, 6)
    plin = heatmap(mf.components; c = :viridis, xticks = xtk, yticks = ytk,
                   title = "linear power — the skirts are invisible", titlefontsize = 10,
                   xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                   left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    pdb = heatmap(db.components; c = :viridis, xticks = xtk, yticks = ytk, clims = (-40, 0),
                  title = "pow2db → still a $(typeof(db).name.name), axes intact: $(db.fcs == mf.fcs)",
                  titlefontsize = 10,
                  xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                  left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    keep(plot(plin, pdb; layout = (1,2), size = (1100, 460)), "step4c-element-ops")
end
display(step4c)


## ---------------------------------------------------------------------------------------
## Step 5 — power_envelope
## ---------------------------------------------------------------------------------------

const vfs = 3125.0
vbank(; kwargs...) = gammatone_filterbank(vfs; fmin = 100, fmax = 400, base_frq = 200,
                                          filters_per_erb = 2, bw_scale = 2, kwargs...)
vam(A, m, fc, fm, dur) = (t = (0:round(Int, dur*vfs)-1)./vfs;
                          A.*(1 .+ m.*cos.(2π*fm.*t)).*cos.(2π*fc.*t))

#5a. Does the envelope actually track the waveform? The gammatone channel is analytic, so the
#envelope is |y| directly - no rectify-and-smooth, no Hilbert transform. With the gain-2
#normalisation it comes out in the units of the signal: an amplitude envelope reads A(1+m cos)
#and a power envelope its square. Both are drawn against the truth, so a gain error of any kind
#shows up as a gap rather than having to be inferred from a number.
step5a = let
    A, m, fc, fm = 1.0, 0.5, 200.0, 7.0
    s = signal(vam(A, m, fc, fm, 20.0), vfs)
    fb = vbank()
    k = frqindex(fc, fb.fcs)

    ta = power_envelope(s, fb; envelope = :amplitude, fs_env = vfs)   #undecimated, to overlay on the waveform
    tp = power_envelope(s, fb; envelope = :power, fs_env = vfs)
    win = 1:round(Int, 2.5*vfs/fm)                                    #two and a half modulation periods
    t = ta.time[win]

    #The time axis is the OUTPUT sample's time, as in gammatone_analysis, so the analysis delay
    #is still in it: an event shows up later by the alignment target nd, plus the difference
    #between the filter's group and envelope delays. Shift the truth by that to compare shapes.
    lagsamp = gammatone_delay(fb; delay = :auto).nd +
              group_delay(fb.filters[k]) - envelope_delay(fb.filters[k])
    truth = A.*(1 .+ m.*cos.(2π*fm.*(t .- lagsamp/vfs)))

    p1 = plot(t, real(vam(A, m, fc, fm, 20.0)[round.(Int, t.*vfs) .+ 1]);
              lc = :grey80, lw = 0.6, label = "signal", legend = :topright,
              ylabel = "amplitude", title = "amplitude envelope rides the peaks (gain-2: reads A, not A/2 or 2A)",
              titlefontsize = 9)
    plot!(p1, t, ta.components[k, win]; lc = :crimson, lw = 2, label = "|y| envelope")
    plot!(p1, t, truth; lc = :black, ls = :dash, lw = 1.5, label = "A(1 + m cos), delayed $(round(lagsamp, digits = 1)) samples")

    p2 = plot(t, tp.components[k, win]; lc = :steelblue, lw = 2, label = "|y|² envelope",
              legend = :topright, xlabel = "Time (s)", ylabel = "power",
              title = "power envelope = the square: mean A²(1+m²/2), swings (1∓m)²",
              titlefontsize = 9)
    plot!(p2, t, truth.^2; lc = :black, ls = :dash, lw = 1.5, label = "A²(1 + m cos)², delayed")
    hline!(p2, [A^2*(1 + m^2/2)]; lc = :green, ls = :dot, label = "mean A²(1+m²/2)")
    hline!(p2, [A^2*(1-m)^2, A^2*(1+m)^2]; lc = :orange, ls = :dot, label = "A²(1∓m)²")
    keep(plot(p1, p2; layout = (2,1), size = (950, 640)), "step5a-envelope-tracks-waveform")
end
display(step5a)

#5b. Decimation is only legitimate because the resampler lowpasses before it decimates. Left: the
#envelope spectrum before and after decimating to 125 Hz - they lie on top of each other through
#the passband. Right: the trap it avoids. A 100 Hz modulation is above the 62.5 Hz envelope
#Nyquist, so a naive take-every-Mth-sample folds it to |125-100| = 25 Hz and invents a line that
#is not in the signal; the anti-aliased path leaves nothing there. The dashed marker is where the
#alias would sit.
step5b = let
    fb1k = gammatone_filterbank(vfs, [1000.0]; bw = [400.0])
    fs_env = vfs/25

    slow = signal(vam(1.0, 0.5, 1000.0, 25.0, 40.0), vfs)
    e_full = power_envelope(comp(), slow, fb1k; fs_env = vfs, trim = false)[1,:]
    e_dec  = power_envelope(comp(), slow, fb1k; fs_env)[1,:]
    Pf, ff = welch_psd(comp(), e_full, vfs; frame_dur = 4.0)
    Pd, fd = welch_psd(comp(), e_dec, fs_env; frame_dur = 4.0)
    pl = plot(ff, pow2db.(Pf); lc = :grey60, lw = 1, label = "envelope at $(vfs) Hz",
              xlims = (0, 62.5), xlabel = "Modulation Rate (Hz)", ylabel = "Power (dB)",
              title = "decimation preserves the passband", titlefontsize = 9, legend = :topright,
              left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    plot!(pl, fd, pow2db.(Pd); lc = :crimson, lw = 1.6, ls = :dash, label = "decimated to $(fs_env) Hz")

    fast = signal(vam(1.0, 0.5, 1000.0, 100.0, 40.0), vfs)
    e_naive = power_envelope(comp(), fast, fb1k; fs_env = vfs, trim = false)[1, 1:25:end]  #no AA filter
    e_aa    = power_envelope(comp(), fast, fb1k; fs_env)[1,:]                              #resampled
    Pn, fn = welch_psd(comp(), e_naive, fs_env; frame_dur = 4.0)
    Pa, fa = welch_psd(comp(), e_aa, fs_env; frame_dur = 4.0)
    pr = plot(fn, pow2db.(Pn); lc = :orange, lw = 1.6, label = "naive: every 25th sample",
              xlims = (0, 62.5), xlabel = "Modulation Rate (Hz)", ylabel = "Power (dB)",
              title = "100 Hz modulation: the alias that anti-aliasing removes",
              titlefontsize = 9, legend = :topright,
              left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    plot!(pr, fa, pow2db.(Pa); lc = :seagreen, lw = 1.6, label = "power_envelope (anti-aliased)")
    vline!(pr, [25.0]; lc = :red, ls = :dash, label = "|125 − 100| = 25 Hz")
    keep(plot(pl, pr; layout = (1,2), size = (1150, 480)), "step5b-decimation-and-aliasing")
end
display(step5b)

#5c. What alignment is for. A click through the bank, drawn as channels x time. Uncompensated,
#each channel responds at its own delay and the click smears into a diagonal sweep - low
#channels are slowest. Compensated, every channel peaks at one instant and the click is a single
#vertical ridge. Per-channel modulation spectra do not care (|DFT|^2 is shift invariant), but
#anything comparing channels does, which is why align defaults to :auto.
step5c = let
    fb = gammatone_filterbank(vfs; fmin = 80, fmax = 1200, base_frq = 308, filters_per_erb = 2, bw_scale = 2)
    click = zeros(2000); click[400] = 1.0
    s = signal(click, vfs)

    E0 = power_envelope(comp(), s, fb; fs_env = vfs, trim = false, align = nothing)
    E1 = power_envelope(comp(), s, fb; fs_env = vfs, trim = false)
    lead = compensation_lead(gammatone_delay(fb; delay = :auto))
    win = 350:600
    ytk = TimeFrequencyAnalysis.generate_ticks(fb.fcs, 6)

    p0 = heatmap(amp2db.(E0[:, win] ./ maximum(E0)); c = :viridis, clims = (-40, 0), colorbar = false,
                 yticks = ytk, xguide = "sample", yguide = "Carrier Frequency (Hz)",
                 title = "align = nothing: each channel at its own delay", titlefontsize = 9,
                 left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    p1 = heatmap(amp2db.(E1[:, win] ./ maximum(E1)); c = :viridis, clims = (-40, 0), colorbar = false,
                 yticks = ytk, xguide = "sample", yguide = "Carrier Frequency (Hz)",
                 title = "align = :auto: one instant, lead of $lead already dropped", titlefontsize = 9,
                 left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)

    sp0 = [argmax(view(E0, k, :)) for k in eachindex(fb.fcs)]
    sp1 = [argmax(view(E1, k, :)) for k in eachindex(fb.fcs)]
    pp = plot(sp0, 1:length(fb.fcs); lc = :orange, lw = 2, label = "align = nothing",
              xlabel = "sample of peak", ylabel = "channel", legend = :bottomright,
              title = "peak instant per channel: spread $(maximum(sp0)-minimum(sp0)) → $(maximum(sp1)-minimum(sp1)) samples",
              titlefontsize = 9)
    plot!(pp, sp1, 1:length(fb.fcs); lc = :seagreen, lw = 2, label = "align = :auto")
    keep(plot(p0, p1, pp; layout = @layout([a b; c]), size = (1100, 760)), "step5c-alignment")
end
display(step5c)

#5d. The rule that decides how much of the modulation axis is real: a channel of bandwidth B
#passes both sidebands fc +/- f_m only while 2*f_m is inside B, so a modulation at f_m simply is
#not visible in a channel much narrower than 2*f_m. No citation for this was found in the
#literature, so it is measured here rather than asserted - and it is what the bandwidth floor in
#SONAR's carrier grid rests on.
step5d = let
    fc, m = 500.0, 0.5
    ratios = exp2.(range(-1.5, 3.5; length = 21))         #B/f_m from 0.35 to about 11
    p = plot(; xscale = :log2, xlabel = "channel bandwidth B / modulation rate fₘ",
             ylabel = "recovered modulation depth", legend = :bottomright, size = (950, 520),
             title = "a modulation at fₘ needs a channel of bandwidth B ≳ 2fₘ — and the curve depends only on B/fₘ",
             titlefontsize = 10, left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    for (fm, col) in ((5.0, :navy), (10.0, :seagreen), (20.0, :crimson))
        s = signal(vam(1.0, m, fc, fm, 20.0), vfs)
        d = Float64[]
        for r in ratios
            fb = gammatone_filterbank(vfs, [fc]; bw = [r*fm])
            e = power_envelope(comp(), s, fb; envelope = :amplitude, max_mod_frq = 3fm)[1,:]
            push!(d, (maximum(e) - minimum(e))/(maximum(e) + minimum(e)))
        end
        plot!(p, ratios, d./m; lc = col, lw = 2, marker = :circle, ms = 2.5, label = "fₘ = $fm Hz")
    end
    vline!(p, [2.0]; lc = :red, ls = :dash, lw = 2, label = "B = 2fₘ")
    hline!(p, [1.0]; lc = :grey, ls = :dot, label = "true depth")
    keep(p, "step5d-bandwidth-rule")
end
display(step5d)


## ---------------------------------------------------------------------------------------
## Step 6 — modulation_spectrum
## ---------------------------------------------------------------------------------------

#6a. The headline figure, and the whole point of a two-dimensional representation. Two carriers,
#150 Hz and 600 Hz, modulated at different rates, 3 Hz and 11 Hz. A modulation spectrum alone
#would show two rates but not say which carrier each belonged to; a spectrogram would show two
#carriers but not their modulation. The joint map shows two isolated cells and - the part that
#matters - NOTHING at the two off-diagonal cells that would pair each carrier with the other's
#rate. dc = :remove here, so the carrier axis reads modulation POWER (see 6c).
step6a = let
    fb = gammatone_filterbank(vfs; fmin = 100, fmax = 800, base_frq = 200, filters_per_erb = 3, bw_scale = 2)
    t = (0:round(Int, 60vfs)-1)./vfs
    x = (1 .+ 0.5cos.(2π*3 .*t)).*cos.(2π*150 .*t) .+ (1 .+ 0.5cos.(2π*11 .*t)).*cos.(2π*600 .*t)
    ms = modulation_spectrum(signal(x, vfs), fb; dc = :remove, frame_dur = 10.0)

    jmax = frqindex(20.0, ms.mod_frqs)
    db = pow2db.(ms.components[:, 1:jmax])
    p = heatmap(db .- maximum(db); c = :viridis, clims = (-50, 0),
                xticks = TimeFrequencyAnalysis.generate_ticks(ms.mod_frqs[1:jmax], 8),
                yticks = TimeFrequencyAnalysis.generate_ticks(ms.fcs, 8),
                xguide = "Modulation Rate (Hz)", yguide = "Carrier Frequency (Hz)",
                title = "150 Hz @ 3 Hz and 600 Hz @ 11 Hz — two cells, and the off-diagonals are empty",
                titlefontsize = 10, size = (1000, 520),
                left_margin = 8Plots.mm, bottom_margin = 6Plots.mm)
    for (fc, fm, lab) in ((150.0, 3.0, "→ here"), (600.0, 11.0, "→ here"),
                          (150.0, 11.0, "✗ empty"), (600.0, 3.0, "✗ empty"))
        scatter!(p, [frqindex(fm, ms.mod_frqs[1:jmax])], [frqindex(fc, ms.fcs)];
                 mc = :white, ms = 6, markershape = :circle, markerstrokecolor = :red,
                 markerstrokewidth = 2, alpha = 0.0, label = "")
        annotate!(p, frqindex(fm, ms.mod_frqs[1:jmax]) + 4, frqindex(fc, ms.fcs),
                  text(lab, 8, :left, :white))
    end
    #The third, weaker line on the 150 Hz carrier is the 2fm harmonic of its own 3 Hz modulation
    #(see 6b), not a cross-term: it sits on the RIGHT carrier, at twice the right rate.
    annotate!(p, frqindex(6.0, ms.mod_frqs[1:jmax]) + 4, frqindex(260.0, ms.fcs),
              text("2fₘ harmonic of the\n3 Hz modulation", 7, :left, :white))
    keep(p, "step6a-two-carriers-two-rates")
end
display(step6a)

#6b. Harmonic structure, and the level check. The power envelope of A(1+m cos)cos(2*pi*fc*t) is
#A^2(1 + m cos)^2 = A^2(1 + m^2/2) + 2mA^2 cos + (m^2/2)A^2 cos(2*), so with scaling = :spectrum
#the fm bin reads 2m^2A^4 and the 2fm bin sits (m/4)^2 below it - 18 dB for m = 0.5, and further
#down as m falls. Both predictions are drawn, so the levels are checked rather than admired.
step6b = let
    fb = vbank()
    fc, fm = 200.0, 7.0
    p = plot(; xlabel = "Modulation Rate (Hz)", ylabel = "Power (dB)", xlims = (0, 25),
             legend = :topright, size = (950, 540), titlefontsize = 10,
             title = "one carrier row: the fₘ line at 2m²A⁴ and its 2fₘ harmonic (m/4)² below",
             left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    for (m, col) in zip((0.1, 0.25, 0.5, 1.0), (:navy, :seagreen, :orange, :crimson))
        s = signal(vam(1.0, m, fc, fm, 60.0), vfs)
        ms = modulation_spectrum(s, fb; dc = :remove, frame_dur = 10.0)
        k = frqindex(fc, ms.fcs)
        plot!(p, ms.mod_frqs, pow2db.(ms.components[k,:]); lc = col, lw = 1.5, label = "m = $m")
        scatter!(p, [fm], [pow2db(2*m^2)]; mc = col, ms = 5, markerstrokewidth = 0, label = "")
        scatter!(p, [2fm], [pow2db(2*m^2*(m/4)^2)]; mc = col, ms = 5, markershape = :diamond,
                 markerstrokewidth = 0, label = "")
    end
    annotate!(p, 13.5, pow2db(2*1.0^2) - 3, text("● theory at fₘ = 2m²A⁴\n◆ theory at 2fₘ", 8, :left, :black))
    keep(p, "step6b-harmonic-structure")
end
display(step6b)

#6c. What dc does to the carrier axis - the thing most likely to be misread. :remove gives
#absolute modulation power, so the carrier axis points at the carrier. :index divides each frame
#by its own mean and gives modulation DEPTH, which for a uniformly modulated signal is the same
#in every channel, including channels that barely hear it. The index map is therefore almost flat
#down the carrier axis: that is correct, not a bug, but it means the largest entry of an :index
#map does not mark the carrier. Both agree on the rate.
step6c = let
    fb = vbank()
    fc, fm, m = 200.0, 7.0, 0.5
    s = signal(vam(1.0, m, fc, fm, 60.0), vfs)
    mr = modulation_spectrum(s, fb; dc = :remove, frame_dur = 10.0)
    mi = modulation_spectrum(s, fb; dc = :index,  frame_dur = 10.0)
    j = frqindex(fm, mr.mod_frqs)
    jmax = frqindex(25.0, mr.mod_frqs)
    ytk = TimeFrequencyAnalysis.generate_ticks(mr.fcs, 6)
    xtk = TimeFrequencyAnalysis.generate_ticks(mr.mod_frqs[1:jmax], 6)

    pr = heatmap(pow2db.(mr.components[:,1:jmax]); c = :viridis, clims = (-60, 0), colorbar = false,
                 xticks = xtk, yticks = ytk, xguide = "Modulation Rate (Hz)",
                 yguide = "Carrier Frequency (Hz)", title = "dc = :remove — power, peaked on the carrier",
                 titlefontsize = 9, left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    pi_ = heatmap(pow2db.(mi.components[:,1:jmax]); c = :viridis, clims = (-60, 0), colorbar = false,
                  xticks = xtk, yticks = ytk, xguide = "Modulation Rate (Hz)",
                  yguide = "Carrier Frequency (Hz)", title = "dc = :index — depth, flat down the carrier axis",
                  titlefontsize = 9, left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)

    pc = plot(pow2db.(mr.components[:,j]), 1:length(mr.fcs); lc = :crimson, lw = 2, label = ":remove (power)",
              yticks = ytk, xlabel = "level at fₘ (dB)", ylabel = "Carrier Frequency (Hz)",
              legend = :bottomright, titlefontsize = 9,
              title = "the fₘ column: power spans $(round(Int, pow2db(maximum(mr.components[:,j])/minimum(mr.components[:,j])))) dB, index only $(round(pow2db(maximum(mi.components[:,j])/minimum(mi.components[:,j])), digits=1)) dB")
    plot!(pc, pow2db.(mi.components[:,j]), 1:length(mi.fcs); lc = :seagreen, lw = 2, label = ":index (depth)")
    keep(plot(pr, pi_, pc; layout = @layout([a b; c]), size = (1100, 780)), "step6c-dc-and-the-carrier-axis")
end
display(step6c)

#6d. Cross-check against a spectrogram on the same data, which is the check that will be run on
#real beams: AM puts sidebands at fc +/- fm, so the sideband spacing a spectrogram shows must
#equal the modulation rate the modulation spectrum reports. Two independent measurements of the
#same physical quantity - if they disagree, the chain is wrong somewhere.
step6d = let
    fb = vbank()
    fc, fm, m = 200.0, 7.0, 0.6
    s = signal(vam(1.0, m, fc, fm, 60.0) .+ 0.05.*randn(round(Int, 60vfs)), vfs)

    sp = welch_psd(s; frame_dur = 8.0, scaling = :spectrum)
    lo, hi = frqindex(170.0, sp.frqs), frqindex(230.0, sp.frqs)
    ps = plot(sp.frqs[lo:hi], pow2db.(sp.components[lo:hi]); lc = :steelblue, lw = 1.4, label = "",
              xlabel = "Frequency (Hz)", ylabel = "Power (dB)", titlefontsize = 9,
              title = "spectrum: carrier with sidebands", left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    vline!(ps, [fc - fm, fc + fm]; lc = :red, ls = :dash, label = "fc ± fₘ")

    ms = modulation_spectrum(s, fb; dc = :remove, frame_dur = 10.0)
    k = frqindex(fc, ms.fcs)
    jm = frqindex(25.0, ms.mod_frqs)
    peak = ms.mod_frqs[argmax(ms.components[k, 2:jm]) + 1]
    pm = plot(ms.mod_frqs[1:jm], pow2db.(ms.components[k,1:jm]); lc = :seagreen, lw = 1.4, label = "",
              xlabel = "Modulation Rate (Hz)", ylabel = "Power (dB)", titlefontsize = 9,
              title = "modulation spectrum reports fₘ = $(round(peak, digits = 2)) Hz; sidebands sit ±$(round(peak, digits = 2)) Hz",
              left_margin = 8Plots.mm, bottom_margin = 5Plots.mm)
    vline!(pm, [fm]; lc = :red, ls = :dash, label = "true fₘ")
    keep(plot(ps, pm; layout = (2,1), size = (950, 640)), "step6d-spectrogram-cross-check")
end
display(step6d)


## ---------------------------------------------------------------------------------------
## Step 7 — the modfreq plot recipe
## ---------------------------------------------------------------------------------------

#7a. The recipe's real job: plot(ms) straight from modulation_spectrum, with nothing set by hand.
#Everything the earlier sections built manually - index-space heatmap so the 116 ERB-spaced
#carriers stay evenly spaced on screen, tick labels carrying the true frequencies, both axes
#named, the title from the object - now comes from the recipe. The bank here is the SONAR default
#(ERB grid over 77-1232 Hz at 8 filters per ERB), so this is what a beam will look like.
step7a = let
    fb = gammatone_filterbank(vfs; fmin = 77, fmax = 1232, base_frq = 308,
                              filters_per_erb = 8, bw_scale = 2)
    t = (0:round(Int, 120vfs)-1)./vfs
    #A shaft line at 4.7 Hz with harmonics on a 300 Hz carrier, a blade line at 23.5 Hz on 800 Hz
    shaft = (1 .+ 0.4cos.(2π*4.7 .*t) .+ 0.2cos.(2π*9.4 .*t) .+ 0.1cos.(2π*14.1 .*t)).*cos.(2π*300 .*t)
    blade = (1 .+ 0.5cos.(2π*23.5 .*t)).*cos.(2π*800 .*t)
    s = signal(shaft .+ blade .+ 0.05.*randn(length(t)), vfs)

    ms = modulation_spectrum(s, fb; dc = :remove, frame_dur = 10.0)
    println("step7a: $(length(ms.fcs)) carriers x $(length(ms.mod_frqs)) modulation bins, " *
            "fs_env = $(ms.fs_env) Hz, $(ms.nframes) frames averaged")

    #Everything below is the recipe; only the colour limits and the modulation range are ours
    jmax = frqindex(40.0, ms.mod_frqs)
    view = modfreq(ms.components[:, 1:jmax], ms.fcs, ms.mod_frqs[1:jmax], ms.bws,
                   ms.fs_env, ms.frame_dur, ms.nframes, ms.title)
    p = plot(pow2db(view); clims = (-70, -10), size = (1000, 560),   #colormap: the recipe's default
             left_margin = 10Plots.mm, bottom_margin = 6Plots.mm)
    keep(p, "step7a-recipe-straight-from-modulation-spectrum")
end
display(step7a)

#7b. The knobs the recipe exposes, and that a dB view stays a modfreq. pow2db returns a modfreq
#with its axes intact, so it plots through the same recipe; nxticks and nyticks are independent,
#which matters because the two axes have very different natural densities - 116 carriers against
#1563 modulation bins.
step7b = let
    fb = gammatone_filterbank(vfs; fmin = 77, fmax = 1232, base_frq = 308,
                              filters_per_erb = 4, bw_scale = 2)
    t = (0:round(Int, 60vfs)-1)./vfs
    s = signal((1 .+ 0.5cos.(2π*7 .*t)).*cos.(2π*400 .*t) .+ 0.05.*randn(length(t)), vfs)
    ms = modulation_spectrum(s, fb; dc = :remove, frame_dur = 10.0)
    jmax = frqindex(25.0, ms.mod_frqs)
    view = modfreq(ms.components[:, 1:jmax], ms.fcs, ms.mod_frqs[1:jmax], ms.bws,
                   ms.fs_env, ms.frame_dur, ms.nframes, ms.title)

    plin = plot(view; colorbar = false, titlefontsize = 9,
                title = "plot(ms) — linear power, default :ice colormap")
    pdb  = plot(pow2db(view); clims = (-80, -10), colorbar = false, titlefontsize = 9,
                nxticks = 5, nyticks = 4, title = "plot(pow2db(ms)); nxticks = 5, nyticks = 4")
    keep(plot(plin, pdb; layout = (1,2), size = (1150, 480),
              left_margin = 8Plots.mm, bottom_margin = 6Plots.mm),
         "step7b-recipe-knobs-and-db")
end
display(step7b)


## ---------------------------------------------------------------------------------------
## Step 7 (continued) — the same contact in noise, and how the colormap changes what you see
## ---------------------------------------------------------------------------------------

#A shaft line at 4.7 Hz with its harmonic comb on a 300 Hz carrier, plus a blade line at 23.5 Hz
#on an 800 Hz carrier: the contact used in 7a, reused here so the only thing changing is noise.
function noisy_contact(snr_db, dur = 120.0)
    t = (0:round(Int, dur*vfs)-1)./vfs
    shaft = (1 .+ 0.4cos.(2π*4.7 .*t) .+ 0.2cos.(2π*9.4 .*t) .+ 0.1cos.(2π*14.1 .*t)).*cos.(2π*300 .*t)
    blade = (1 .+ 0.5cos.(2π*23.5 .*t)).*cos.(2π*800 .*t)
    clean = shaft .+ blade
    pc = sum(abs2, clean)/length(clean)
    σ = sqrt(pc/exp10(snr_db/10))
    return signal(clean .+ σ.*randn(length(t)), vfs)
end

sonar_bank() = gammatone_filterbank(vfs; fmin = 77, fmax = 1232, base_frq = 308,
                                    filters_per_erb = 8, bw_scale = 2)

#Trim a modfreq to a modulation range, keeping it a modfreq so the recipe still draws it
function modview(ms, fmax)
    j = frqindex(fmax, ms.mod_frqs)
    return modfreq(ms.components[:, 1:j], ms.fcs, ms.mod_frqs[1:j], ms.bws,
                   ms.fs_env, ms.frame_dur, ms.nframes, ms.title)
end

#Express a modfreq in dB above its own median floor. The absolute level of a modulation spectrum
#rises with the noise, so a fixed colour scale saturates as soon as the noise does; what a reader
#actually wants to know is how far a line stands above the floor around it, which is what a sonar
#display shows and what a detector thresholds.
function floor_db(v)
    db = pow2db.(v.components)
    return modfreq(db .- median(db), v.fcs, v.mod_frqs, v.bws, v.fs_env, v.frame_dur, v.nframes, v.title)
end

#7c. The same contact at falling broadband SNR, each panel scaled to its own noise floor. Two
#things to notice. First, the lines survive far below 0 dB BROADBAND SNR: a gammatone channel is
#narrow, so almost none of the band's noise reaches any one channel, and the shaft line still
#stands 18 dB clear at -10 dB SNR. Second, noise does not blur or move the lines - broadband noise
#through a channel gives a fluctuating envelope, which is a raised pedestal across ALL modulation
#rates - so the lines stay exactly where they are and are simply overtaken. The weak third shaft
#harmonic goes first, the blade line last.
step7c = let
    fb = sonar_bank()
    panels = Any[]
    for snr in (20, 0, -10, -20)
        v = floor_db(modview(modulation_spectrum(noisy_contact(snr), fb; dc = :remove, frame_dur = 10.0), 40.0))
        k, jf = frqindex(300.0, v.fcs), frqindex(4.7, v.mod_frqs)
        kb, jb = frqindex(800.0, v.fcs), frqindex(23.5, v.mod_frqs)
        push!(panels, plot(v; clims = (0, 40), colorbar = false,     #no c = : the recipe's :ice default
                           nxticks = 5, nyticks = 4, titlefontsize = 9,
                           title = "SNR $(snr) dB — shaft $(round(v.components[k,jf], digits = 1)) dB, blade $(round(v.components[kb,jb], digits = 1)) dB over floor"))
    end
    keep(plot(panels...; layout = (2,2), size = (1200, 800),
              left_margin = 8Plots.mm, bottom_margin = 5Plots.mm), "step7c-contact-in-noise")
end
display(step7c)

#7d. The same noisy data (SNR -10 dB), the same colour limits, six colormaps. Only the colormap
#changes, so every difference is the colormap talking. This is SEQUENTIAL data - a magnitude - and
#the rule for magnitude is one hue running light to dark. The bottom row breaks that rule, and the
#cost is that lightness stops encoding magnitude: turbo and jet put their brightest colour in the
#MIDDLE of the range, so the mid-level NOISE PEDESTAL is rendered as brightly as the lines and
#competes with them for attention. In the top row the floor stays near-black and only the contact
#is bright. 7e measures the cause and shows what each map looks like with the colour taken away.
step7d = let
    fb = sonar_bank()
    v = floor_db(modview(modulation_spectrum(noisy_contact(-10), fb; dc = :remove, frame_dur = 10.0), 40.0))
    panels = [plot(v; c = cm, clims = (0, 20), colorbar = false, nxticks = 4, nyticks = 4,
                   titlefontsize = 10,
                   title = "$cm" * (cm in (:turbo, :jet) ? "  ✗ rainbow" : ""))
              for cm in (:viridis, :magma, :ice, :greys, :turbo, :jet)]
    keep(plot(panels...; layout = (2,3), size = (1300, 740),
              left_margin = 7Plots.mm, bottom_margin = 5Plots.mm), "step7d-colormap-comparison")
end
display(step7d)

#7e. Why the bottom row of 7d misleads, measured rather than asserted. A colormap encodes
#magnitude honestly only if perceived lightness rises monotonically across it: equal steps in the
#data must look like equal steps on the page. Left: L* (CIE lightness) along each map. The
#perceptual maps climb steadily from dark to light. jet and turbo do not - they start dark, peak
#bright in the MIDDLE, and end dark again, so their two ends are indistinguishable in lightness
#and a mid-range value outshines the maximum. Right: the same data drawn as a greyscale of each
#map's lightness alone, which is what a colourblind reader, a photocopy or a bad projector sees.
step7e = let
    C = Plots.Colors
    Lstar(cm, n = 128) = [C.convert(C.Lab, C.RGB(get(cgrad(cm), t))).l for t in range(0, 1; length = n)]

    pl = plot(; xlabel = "position along the colormap", ylabel = "perceived lightness L*",
              legend = :bottomright, titlefontsize = 10, ylims = (-2, 102),
              title = "honest magnitude needs L* to rise monotonically")
    for (cm, col, st) in ((:viridis, :seagreen, :solid), (:magma, :purple, :solid),
                          (:ice, :steelblue, :solid), (:greys, :black, :solid),
                          (:turbo, :orange, :dash), (:jet, :crimson, :dash))
        L = Lstar(cm)
        mono = all(diff(L) .> -0.5)
        plot!(pl, range(0, 1; length = length(L)), L; lc = col, ls = st, lw = 2,
              label = "$cm$(mono ? "" : "  ✗ not monotonic")")
    end

    fb = sonar_bank()
    v = floor_db(modview(modulation_spectrum(noisy_contact(-10), fb; dc = :remove, frame_dur = 10.0), 40.0))
    #Redraw each map using ONLY its lightness: what survives if colour is taken away
    greyof(cm) = cgrad([C.RGB(l/100, l/100, l/100) for l in Lstar(cm, 64)])
    gp = [plot(v; c = greyof(cm), clims = (0, 20), colorbar = false,
               nxticks = 3, nyticks = 3, titlefontsize = 9, title = "$cm — lightness only")
          for cm in (:viridis, :magma, :turbo, :jet)]
    keep(plot(pl, plot(gp...; layout = (2,2)); layout = (1,2), size = (1350, 560),
              left_margin = 7Plots.mm, bottom_margin = 5Plots.mm), "step7e-why-lightness-matters")
end
display(step7e)
