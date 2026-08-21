```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Filter Responses

A rational digital filter is represented by its numerator/denominator coefficients in a
[`filter_coefs`](@ref). Signals are filtered with `filt` — DSP.jl's filtering function,
re-exported with methods for the package's filter types (`filter_coefs` here, and the
[gammatone types](gammatone.md)) so that generic code can filter through one verb. Each
`filt`/`filt!` pair also has a **stateful form** for frame-by-frame (streaming) work: it
takes a filter state — created at rest by [`filter_state`](@ref) — as an extra argument,
starts from it instead of from zero initial conditions, and leaves the end-of-block state
in it for the next call, so consecutive chunks filter as if the signal had been processed
whole. For a `filter_coefs` the state follows DSP.jl's `DF2TFilter` convention exactly and
the two are interchangeable mid-stream. The
complex and magnitude frequency responses are evaluated on arbitrary frequency grids with
[`filter_resp`](@ref) and [`filter_magresp`](@ref) — for example, to compute the spectrum
of an all-pole (AR) model over a chosen set of analysis frequencies — and the impulse
response is measured in the time domain with [`filter_impresp`](@ref), which works for
any filter with a `filt` method. The pointwise evaluation and frequency-conversion
helpers behind these are documented under [Internals](internals.md).

```@docs
filter_coefs
filt(::filter_coefs, ::AbstractVector)
filt!(::AbstractVector, ::filter_coefs, ::AbstractVector)
filter_state
filt(::filter_coefs, ::AbstractVector, ::AbstractVector)
filt!(::AbstractVector, ::filter_coefs, ::AbstractVector, ::AbstractVector)
filter_resp
filter_magresp
filter_impresp
```
