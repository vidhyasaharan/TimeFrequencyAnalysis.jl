```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Correlations

Autocorrelation and cross-correlation of discrete-time sequences, plus short-time
correlation functions evaluated directly on a [`signal`](@ref).

## Cross-correlation sequences

[`xcorr`](@ref) and [`xcorr!`](@ref) compute the cross correlation between two sequences —
typically used to search along a longer sequence for the best match to a shorter one.
Implemented as a *sliding dot product* with zero padding.

```@docs
xcorr
xcorr!
```

## Autocorrelation sequences

[`acorr`](@ref) and [`acorr!`](@ref) compute the autocorrelation sequence up to the
desired lag. The value at each lag is normalised by the size of the correlation window,
which becomes smaller at the edges.

```@docs
acorr
acorr!
```

## Short-time correlation functions

[`acf`](@ref) evaluates the autocorrelation function of a signal at a specified position,
at one or more lags, over a window of specified size. The normalised variant,
[`nacf`](@ref), mean-subtracts both windows and scales the product by their norms so
values lie in `[-1, 1]`; it is the core of correlation-based pitch estimators.

```@docs
acf
nacf
```
