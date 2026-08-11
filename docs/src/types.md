```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Data Structures

Analyses in this package take either a signal or frames of a signal as input. Two structs
serve this purpose:

- [`signal`](@ref) — signal samples together with the sampling rate
- [`framed_signal`](@ref) — a signal plus the geometry for splitting it into frames

Spectral analyses store their output in one of two container structs, which keep the
provenance of the result (the analysed signal, the frequency and time axes, and a title):

- [`spectrum`](@ref) — spectral components over a frequency axis
- [`timefreq`](@ref) — spectro-temporal components over time and frequency axes

Every struct is parametric in `T<:AbstractFloat`, the precision of the signal; constructors
convert all numeric fields to that precision.

## `signal`

```@docs
signal
```

## `framed_signal`

```@docs
framed_signal
```

## `spectrum`

```@docs
spectrum
```

## `timefreq`

```@docs
timefreq
```

## The `comp()` convention

```@docs
comp
```
