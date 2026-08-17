```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Plotting

Plot recipes for [`signal`](@ref), [`spectrum`](@ref), [`timefreq`](@ref) and
[`modfreq`](@ref) are provided by a package extension that loads automatically when RecipesBase is in the
session — in practice, whenever Plots is loaded. There is no extra installation step:

```julia
using Plots, TimeFrequencyAnalysis

s = signal(x, fs)
plot(s)                                  #signal against time in seconds
plot(magspec(s))                         #spectrum against its frequency grid
plot(amp2db(specgram(framed_signal(s)))) #spectrogram heatmap in dB
plot(pow2db(modulation_spectrum(s, fb))) #modulation spectrum heatmap in dB
```

- A `signal` plots as a line of amplitude against time in seconds.
- A `spectrum` plots as a line against its frequency grid, with `nticks` (default 8)
  controlling the number of frequency ticks.
- A `timefreq` plots as a heatmap, with `nxticks`/`nyticks` (default 8) controlling the
  tick counts. Heatmaps are drawn in index space with ticks labelled by the true
  time/frequency coordinates: index space renders nonuniform (e.g. logarithmic) frequency
  grids uniformly.
- A `modfreq` plots as a heatmap on the same terms, with carrier frequency on the vertical
  axis and modulation rate on the horizontal one.

No dB or other scaling is applied implicitly — pass [`amp2db`](@ref) etc. explicitly.

Everything the recipes set apart from the series type is a **default**, so any of it can be
overridden in the usual way — `plot(tf; title = "…", xguide = "…", xticks = …, c = …)`. The
series type is forced: a `timefreq` or `modfreq` is a heatmap.

## Colormaps

Both heatmap recipes default to **`:ice`**. That is a deliberate choice, and worth understanding
before overriding it.

A spectrogram, cochleagram or modulation spectrum is a **magnitude** — an ordered quantity with
no meaningful midpoint. The colormap for a magnitude must be *sequential*: one hue running dark
to light, with **perceived lightness rising monotonically** across the map, so that equal steps
in the data look like equal steps on the page. Lightness is what the eye reads as "more"; hue is
decoration.

Measured lightness (CIE `L*`, 0 = black, 100 = white) across the maps most often reached for:

| colormap | `L*` start → end | lightness monotonic? | |
|:---|:---|:---|:---|
| `:ice` | 2 → 98 | yes | **default** — near-black floor, blue through white |
| `:magma`, `:inferno` | 0 → 98 | yes | highest contrast; black floor through red to white |
| `:greys` | 0 → 100 | yes | for print, and the honest baseline |
| `:thermal`, `:haline` | 13–17 → 95 | yes | oceanographic (cmocean), a lighter floor |
| `:viridis`, `:plasma` | 15 → 91–94 | yes | fine, but the floor is a visible purple |
| `:deep` | 98 → 12 | reversed | light-to-dark; use for a white-background figure |
| `:turbo` | 12 → **24** | **NO** | avoid |
| `:jet` | 13 → **25** | **NO** | avoid |

**Why the last two are excluded.** A rainbow map is not monotonic in lightness: `jet` and `turbo`
both start dark, peak bright in the *middle* of their range, and return to dark at the top. Two
consequences follow, and neither is cosmetic. A mid-level value — a noise pedestal, say — is
rendered *more* prominently than the maximum it should be subordinate to, so the eye is drawn to
the wrong thing. And because the two ends sit at nearly the same lightness (12 and 24 for
`turbo`), the strongest and weakest parts of the image are indistinguishable to anyone reading it
in greyscale: a colourblind viewer, a photocopy, a projector with the contrast wound down. The
same maps also invent banding where the data is smooth, because their lightness gradient is
uneven.

`TimeFrequencyAnalysis/scripts/modulation-verification.jl` (sections 7d and 7e) demonstrates all
of this on a real modulation spectrum: the same data under six maps, the lightness curves
measured, and each map redrawn using its lightness alone.

**Choosing one.** For a display where faint structure sits on a noise floor — which is the usual
case — pick a map whose floor is near-black, so the floor recedes and only the signal carries
brightness: `:ice`, `:magma`, `:inferno` or `:greys`. `:viridis` is perfectly sound but bottoms
out at `L* = 15`, so its floor stays a visible purple that competes slightly with weak features.
For a figure on a white page, `:deep` runs the other way and reads better in print.

```julia
plot(pow2db(ms))                #:ice, the default
plot(pow2db(ms); c = :magma)    #or any other sequential map
plot(pow2db(ms); c = :greys)    #for print
```

## Helpers

```@docs
TimeFrequencyAnalysis.generate_ticks
```
