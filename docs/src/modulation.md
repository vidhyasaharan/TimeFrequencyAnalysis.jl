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
using TimeFrequencyAnalysis, Plots

#A bank over the band of interest. bw_scale = 2 widens each channel to twice its ERB, which
#is what lets the sidebands of a modulation through - see the note on bandwidth below.
fb  = gammatone_filterbank(fs; fmin = 100, fmax = 1200, base_frq = 300,
                           filters_per_erb = 4, bw_scale = 2)

env = power_envelope(s, fb)              #timefreq: channels x time, decimated
ms  = modulation_spectrum(s, fb)         #modfreq: carrier x modulation rate
plot(pow2db(ms))                         #heatmap, see the plotting page
```

!!! note "Channel bandwidth limits the modulation rates you can see"
    A channel of bandwidth `B` passes both sidebands of a modulation at `f_m` only while
    `2f_m ≲ B`, so a narrow channel is blind to fast modulation however finely the envelope
    is sampled. Each result carries the per-channel bandwidths in its `bws` field for exactly
    this reason: entry `[k, j]` is only meaningful while `mod_frqs[j] ≲ bws[k]/2`. Widen the
    channels with `bw_scale` (or set `bw` directly) when the modulation rates of interest are
    fast relative to the carrier spacing.

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
