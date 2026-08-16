```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Filter Responses

A rational digital filter is represented by its numerator/denominator coefficients in a
[`filter_coefs`](@ref); its complex and magnitude frequency responses are evaluated on
arbitrary frequency grids with [`filter_resp`](@ref) and [`filter_magresp`](@ref) — for
example, to compute the spectrum of an all-pole (AR) model over a chosen set of analysis
frequencies — and its impulse response is measured in the time domain with
[`filter_impresp`](@ref). The pointwise evaluation and frequency-conversion helpers behind
these are documented under [Internals](internals.md).

```@docs
filter_coefs
filter_resp
filter_magresp
filter_impresp
```
