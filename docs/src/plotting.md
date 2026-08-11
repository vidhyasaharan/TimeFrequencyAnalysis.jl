```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Plotting

Plot recipes for [`signal`](@ref), [`spectrum`](@ref) and [`timefreq`](@ref) are
provided by a package extension that loads automatically when RecipesBase is in the
session — in practice, whenever Plots is loaded. There is no extra installation step:

```julia
using Plots, TimeFrequencyAnalysis

s = signal(x, fs)
plot(s)                                  #signal against time in seconds
plot(magspec(s))                         #spectrum against its frequency grid
plot(amp2db(specgram(framed_signal(s)))) #spectrogram heatmap in dB
```

- A `signal` plots as a line of amplitude against time in seconds.
- A `spectrum` plots as a line against its frequency grid, with `nticks` (default 8)
  controlling the number of frequency ticks.
- A `timefreq` plots as a heatmap, with `nxticks`/`nyticks` (default 8) controlling the
  tick counts. Heatmaps are drawn in index space with ticks labelled by the true
  time/frequency coordinates: index space renders nonuniform (e.g. logarithmic) frequency
  grids uniformly.

No dB or other scaling is applied implicitly — pass [`amp2db`](@ref) etc. explicitly.

## Helpers

```@docs
TimeFrequencyAnalysis.generate_ticks
```
