```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Spectral Analyses

Spectral (frequency axis) and spectro-temporal (time and frequency axes) analyses:

- *Fourier spectrum* — [`dft`](@ref) and [`magspec`](@ref), based on the DFT
- *STFT spectrogram* — [`specgram`](@ref), the short-time Fourier transform magnitude
- *Periodograms* — [`periodogram`](@ref), projections onto complex exponentials at an
  arbitrary frequency grid (for example, log-spaced grids from [`logfreq_array`](@ref))
- *Welch power spectral density* — [`welch_psd`](@ref), frame-averaged power spectrum

Each analysis returns a [`spectrum`](@ref) or [`timefreq`](@ref) container, or the plain
component array when called with [`comp()`](@ref comp) as the first argument.

## Fourier Spectrum

```@docs
dft
magspec
```

## STFT Spectrogram

```@docs
specgram
```

## Periodogram

```@docs
periodogram
```

## Welch Power Spectral Density

```@docs
welch_psd
```

## Working in dB and log scales

`log`, `log10`, `amp2db` and `pow2db` apply elementwise to [`spectrum`](@ref) and
[`timefreq`](@ref) objects and return a new object of the same kind, leaving the signal,
axes and title untouched (`amp2db` and `pow2db` extend the functions of the same name from
DSP.jl):

```julia
tf = specgram(frames)
db = amp2db(tf)      #still a timefreq, with components in dB
lg = log10(magspec(s))
```

```@docs
log(::spectrum)
log10(::spectrum)
amp2db(::spectrum)
pow2db(::spectrum)
```
