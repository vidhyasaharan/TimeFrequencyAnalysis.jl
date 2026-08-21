#Visual verification companion to the stateful (frame-by-frame) filtering step: the
#filt/filt! methods that take a filter_state, carry it across calls and update it in
#place, for gammatone filters, filterbanks and filter_coefs (DF2T). Each section produces
#a plot that shows one property the tests assert numerically. Run with Plots installed
#(Plots is not a package dependency): julia> include("scripts/stateful-filtering-verification.jl")
#Not part of the package or its test suite.

using TimeFrequencyAnalysis
import DSP
using Plots
using Random

## ---------------------------------------------------------------------------------------
## Section 1 — what the carried state does: seamless frames vs per-frame restarts
## ---------------------------------------------------------------------------------------
#A 1 kHz tone filtered by the 1 kHz gammatone channel, processed in 20 ms frames. With the
#state carried between calls the framed envelope lies exactly on the whole-signal one
#(difference at the floating-point floor, shown separately); restarting every frame from
#rest instead re-runs the filter's onset transient at every boundary — the sawtoothed
#envelope is what the stateful form exists to prevent.
section1 = let
    fs = 16000.0
    gf = gammatone_filter(fs, 1000.0)
    t = (0:3199)./fs                      #200 ms
    x = sin.(2π*1000.0.*t)
    L = 320                               #20 ms frames
    frames = [x[(1+(i-1)*L):(i*L)] for i in 1:length(x)÷L]

    ywhole = filt(gf, x)
    w = filter_state(gf)
    ystream = reduce(vcat, [filt(gf, frame, w) for frame in frames])
    yrestart = reduce(vcat, [filt(gf, frame) for frame in frames])

    p1 = plot(1000 .* t, abs.(ywhole); lw = 3, label = "whole signal",
              ylabel = "envelope", title = "20 ms frames through the 1 kHz channel")
    plot!(p1, 1000 .* t, abs.(ystream); lw = 1, label = "framed, state carried")
    plot!(p1, 1000 .* t, abs.(yrestart); label = "framed, restarted from rest")
    vline!(p1, 1000 .* (L:L:length(x)-1)./fs; ls = :dot, lc = :grey, label = "frame boundaries")

    err = abs.(ystream .- ywhole)
    p2 = plot(1000 .* t, err; label = "", xlabel = "Time (ms)",
              ylabel = "|framed - whole|",
              title = "carried state: difference is identically zero (max = $(maximum(err)))")

    plot(p1, p2; layout = (2,1), size = (850, 600))
end
display(section1)

## ---------------------------------------------------------------------------------------
## Section 2 — the filterbank stream: a chirp cochleagram assembled frame by frame
## ---------------------------------------------------------------------------------------
#The linear-chirp cochleagram of the docs page, computed three ways: whole signal, in
#50 ms blocks with per-channel states carried (identical to the last bit), and in 50 ms
#blocks restarted from rest — where every block boundary re-triggers every channel's
#onset, striping the picture. The middle and left panels being indistinguishable is the
#point; the difference panel confirms it is exact.
section2 = let
    fs = 16000.0
    fb = gammatone_filterbank(fs)
    t = (0:Int(fs)-1)./fs
    chirp = sin.(2π.*(100.0.*t .+ 3450.0.*t.^2))  #100 Hz → 7 kHz over one second
    L = 800                                        #50 ms blocks
    blocks = [chirp[(1+(i-1)*L):(i*L)] for i in 1:length(chirp)÷L]

    Ywhole = filt(fb, chirp)
    W = filter_state(fb)
    Ystream = reduce(hcat, [filt(fb, blk, W) for blk in blocks])
    Yrestart = reduce(hcat, [filt(fb, blk) for blk in blocks])

    clims = (-80, 0)
    heat(Y, ttl) = heatmap(t, 1:length(fb.fcs), amp2db.(abs.(Y) .+ eps());
                           clims, title = ttl, ylabel = "channel", colorbar = false)
    p1 = heat(Ywhole, "whole signal")
    p2 = heat(Ystream, "50 ms blocks, state carried")
    p3 = heat(Yrestart, "50 ms blocks, restarted (striped)")
    p4 = heatmap(t, 1:length(fb.fcs), log10.(abs.(Ystream .- Ywhole) .+ 1e-300);
                 title = "log10 |carried - whole| (-300 = exact)", xlabel = "Time (s)",
                 ylabel = "channel", colorbar = true)
    plot(p1, p2, p3, p4; layout = (2,2), size = (1100, 750))
end
display(section2)

## ---------------------------------------------------------------------------------------
## Section 3 — what the state is, and the filter_coefs/DSP hand-off
## ---------------------------------------------------------------------------------------
#Left: after each block of an impulse response run, the carried gammatone state w[j] must
#sit exactly on stage j's closed-form (unnormalised) impulse response
#binomial(n+j-1, j-1)ãⁿ — the block-end markers land on the analytic curves.
#Right: a filter_coefs stream started with the package's stateful filt and handed
#mid-signal — same coefficients, same state vector — to DSP.jl's DF2TFilter continues
#without a seam: the spliced output lies on the whole-signal one.
section3 = let
    fs = 16000.0
    gf = gammatone_filter(fs, 1000.0)
    N = 240
    e = zeros(N); e[1] = 1.0
    L = 40
    w = filter_state(gf)
    ends = Int[]
    states = Vector{ComplexF64}[]
    for i in 1:N÷L
        filt(gf, e[(1+(i-1)*L):(i*L)], w)
        push!(ends, i*L - 1)              #0-based sample index of each block's last sample
        push!(states, copy(w))
    end
    ã = ComplexF64(gf.coef)
    n = 0:N-1
    p1 = plot(; yscale = :log10, xlabel = "sample n", ylabel = "|w_j[n]|",
              title = "carried state vs closed-form stage responses")
    for j in 1:4
        plot!(p1, n, abs.([binomial(k + j - 1, j - 1)*ã^k for k in n]);
              label = "stage $j: binom(n+$(j-1),$(j-1))|ã|ⁿ", lc = j)
        scatter!(p1, ends, abs.([st[j] for st in states]); label = "", mc = j, ms = 5)
    end

    F = filter_coefs([1.0], [1.0, -0.95*cis(2π*500.0/fs)])  #complex one-pole resonator
    x = 0.3 .* randn(Xoshiro(7), 800) .+ sin.(2π*500.0.*(0:799)./fs)
    si = filter_state(F)
    yhead = filt(F, x[1:400], si)
    df = DSP.DF2TFilter(DSP.PolynomialRatio(F.num, F.den), si)  #hand-off: same state vector
    ytail = DSP.filt(df, x[401:800])
    ywhole = filt(F, x)
    p2 = plot((0:799)./fs .* 1000, abs.(ywhole); lw = 3, label = "whole signal (stateless filt)",
              xlabel = "Time (ms)", ylabel = "|y|",
              title = "mid-stream hand-off to DSP.DF2TFilter")
    plot!(p2, (0:399)./fs .* 1000, abs.(yhead); lw = 1, label = "head: TFA stateful filt")
    plot!(p2, (400:799)./fs .* 1000, abs.(ytail); lw = 1, label = "tail: DSP.DF2TFilter, same state")
    vline!(p2, [400/fs*1000]; ls = :dot, lc = :grey, label = "hand-off")

    plot(p1, p2; layout = (1,2), size = (1150, 480))
end
display(section3)

## ---------------------------------------------------------------------------------------
## Section 4 — timings: the register kernel, the general path, and the cost of statefulness
## ---------------------------------------------------------------------------------------
#Printed, not plotted. Answers two design questions on the record: (a) carrying the state
#costs nothing measurable at any frame size (the four states cross memory once per call,
#not per sample), so the stateful form needs no separate implementation; (b) the order-4
#register kernel is ~4× faster per stage than the general stage-sweep path, which is why
#the internal order == 4 split exists — but it serves both entry points, so it stays one
#function with one branch, not two public functions.
section4 = let
    bench(f, reps) = (f(); minimum(@elapsed(f()) for _ in 1:reps))
    fs = 16000.0
    gf4 = gammatone_filter(fs, 1000.0)
    gf5 = gammatone_filter(fs, 1000.0; order = 5)
    x = randn(Xoshiro(1), 160000)
    y = Vector{ComplexF64}(undef, length(x))
    w4 = filter_state(gf4); w5 = filter_state(gf5)
    println("160k samples, order 4:  stateless ", round(bench(() -> filt!(y, gf4, x), 100)*1e3; digits = 3),
            " ms | stateful ", round(bench(() -> filt!(y, gf4, x, w4), 100)*1e3; digits = 3), " ms")
    t5 = bench(() -> filt!(y, gf5, x, w5), 50)
    t4 = bench(() -> filt!(y, gf4, x, w4), 50)
    println("per-stage cost:  order-4 register kernel ", round(t4/4*1e6; digits = 1),
            " µs | general stage sweep ", round(t5/5*1e6; digits = 1), " µs  (×", round(t5/5/(t4/4); digits = 1), ")")
    for L in (160, 512, 2048)
        yf = Vector{ComplexF64}(undef, L); xf = x[1:L]; wf = filter_state(gf4)
        tless = bench(() -> filt!(yf, gf4, xf), 2000)
        tful = bench(() -> filt!(yf, gf4, xf, wf), 2000)
        println("frame L = ", lpad(L, 4), ":  stateless ", round(tless*1e6; digits = 2),
                " µs | stateful ", round(tful*1e6; digits = 2), " µs")
    end
end
