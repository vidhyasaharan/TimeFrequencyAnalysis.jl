#Visual verification companion to the gammatone filterbank plan
#(dev docs/gammatone-filterbank.md §4): one section per implementation step, each producing
#the plots that step's "Verify visually" item calls for. Run with Plots installed (Plots is
#not a package dependency): julia> include("dev docs/gammatone-verification.jl")
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
