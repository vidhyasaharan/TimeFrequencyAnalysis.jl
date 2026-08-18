```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Utilities

Supporting functionality used throughout the package and useful on its own.

## Window functions

```@docs
window
```

## Frequency grids and indexing

Frequency grids define the analysis frequencies for [`periodogram`](@ref) and can be used
to build custom analyses; [`frqindex`](@ref) maps frequencies back to grid indices. The
ERB family implements the auditory (equivalent rectangular bandwidth) frequency scale
behind the [gammatone filterbank](gammatone.md).

```@docs
logfreq_array
linfreq_array
erb
freq2erb
erb2freq
erbfreq_array
frqindex
```

## Reading a number off a signal

Helpers that answer a question about a signal or a spectrum without computing a new
representation of it.

```@docs
peak_frequency
```

## Sampling helpers

```@docs
time2nsamples
resample
```

## Synthetic signal generators

Deterministically shaped random and periodic test signals, used by this package's test
suite and convenient for experiments and examples.

```@docs
white_noise
ar_process
impulse_train
```

## Padding

```@docs
TimeFrequencyAnalysis.zero_pad
TimeFrequencyAnalysis.symmetric_pad
TimeFrequencyAnalysis.unpad_vector
```
