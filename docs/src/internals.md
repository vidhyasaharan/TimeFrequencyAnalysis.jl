```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Internals

Notes for maintainers, plus docstrings of unexported helpers.

## Precision parametrisation

Every container type is parametric in `T<:AbstractFloat`, the element type of the signal it
was built from, and every routine allocates its buffers and outputs in `T`. The rules used
throughout the source:

- literal constants inside generic code are converted with `T(...)` so they do not promote
  `Float32` computations to `Float64`;
- outputs of `Float64`-only library calls (the DSP.jl window shapes, for example) are
  converted once per call with `convert(Vector{T}, ...)`, which is free when `T == Float64`;
- frequency *grids* ([`logfreq_array`](@ref), [`linfreq_array`](@ref)) are always generated
  in `Float64`; the container constructors convert them to `T` on storage.

## Constructor funnelling

`spectrum` and `timefreq` have one central converting constructor that all convenience
constructors call; it converts every field to the signal's precision and then invokes the
default (implicit) constructor. When adding constructors, funnel through it rather than
calling the default constructor directly.

## Unexported helpers

```@docs
findclosest
cexp
cexp_proj_matrix
```

## Filter evaluation helpers

Frequency conversion and pointwise transfer-function evaluation behind
[`filter_resp`](@ref) and [`filter_magresp`](@ref).

```@docs
freq2θ
freq2θ!
freq2z
freq2z!
powers
H
Hmag
```
