```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Modulation Analysis

Modulation analysis describes how the energy in each channel of an auditory filterbank
*fluctuates over time*, rather than how it is distributed across frequency at one instant.
Two functions cover the chain:

- [`power_envelope`](@ref) — the envelope of every channel of a
  [gammatone filterbank](gammatone.md), decimated to a sampling rate suited to the
  modulation rates of interest, returned as a [`timefreq`](@ref)
- [`modulation_spectrum`](@ref) — the power in each of those envelopes at each modulation
  rate, returned as a [`modfreq`](@ref) whose rows are carrier channels and whose columns
  are modulation rates

Both take the signal and a filterbank designed for its sampling rate, and both accept
[`comp()`](@ref comp) as a first argument to get the plain component matrix instead of a
container.

```julia
fb = gammatone_filterbank(fs; nchannels = 32)
env = power_envelope(s, fb)              #timefreq: channels x time
ms  = modulation_spectrum(s, fb)         #modfreq: carrier x modulation rate
plot(pow2db(ms))                         #heatmap, see the plotting page
```

## Envelopes

```@docs
power_envelope
```

## Modulation spectrum

```@docs
modulation_spectrum
```

## Transient allowance

```@docs
startup_samples
```
