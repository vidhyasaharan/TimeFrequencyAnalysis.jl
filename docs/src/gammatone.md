```@meta
CurrentModule = TimeFrequencyAnalysis
```

# Gammatone Filterbank

Auditory-model time-frequency analysis after Hohmann (2002). *Under construction*: the
single-filter design and its inspection helpers ship first; the filterbank, cochleagram
analyses and delay/phase alignment follow (see `dev docs/gammatone-filterbank.md` in the
repository for the full plan).

A [`gammatone_filter`](@ref) is a cascade of `order` (default 4) identical one-pole complex
filters with pole ``ã = λe^{iβ}`` at the centre frequency. Filtering a real signal with it
yields a complex subband signal whose magnitude is the band envelope and whose real part is
the bandpass-filtered signal; the magnitude response peaks at the centre frequency with
gain exactly 2, which is what makes that interpretation hold. The design bandwidth
defaults to the auditory bandwidth [`erb`](@ref)`(fc)` — scaled by `bw_scale`, or replaced
outright via the `bw` keyword (or the explicit-bandwidth constructor form).

The impulse response is available in closed form through [`impulse_response`](@ref), and
[`filter_coefs`](@ref) converts the cascade to its exact rational form so the standard
response machinery applies: [`filter_resp`](@ref) and [`filter_magresp`](@ref) accept a
`gammatone_filter` directly.

Signals are filtered with [`gammatone_filt`](@ref) (or the in-place
[`gammatone_filt!`](@ref)), which applies the cascade as `order` one-pole complex
recursions per sample, computing in the signal's own precision.

A [`gammatone_filterbank`](@ref) collects one filter per channel on an ERB-spaced grid of
centre frequencies (or any explicit grid), and the two analyses build on it:
[`gammatone_analysis`](@ref) returns the complex subband signals as a
[`timefreq`](@ref) — per-sample time resolution, envelope and fine structure per
channel — while [`gammatone_cochleagram`](@ref) returns the envelope form that plots
directly with the existing recipes (`plot(amp2db(gammatone_cochleagram(s)))`). Both follow
the package's [`comp()`](@ref comp) convention for bare-array output.

Gammatone channels are slow at low frequencies and fast at high ones, and two delay
notions matter when interpreting or aligning their outputs: the group delay at the centre
frequency ([`group_delay`](@ref), ≈ `order/(2πb)`) and the slightly earlier peak of the
impulse-response envelope ([`envelope_delay`](@ref), ≈ `(order-1)/(2πb)`), which is what
alignment across channels targets. [`summed_resp`](@ref) shows the bank's summed
frequency response — deeply rippled until the channels are delay- and phase-aligned.

```@docs
gammatone_filter
impulse_response
gammatone_filt
gammatone_filt!
gammatone_filterbank
gammatone_analysis
gammatone_cochleagram
group_delay
envelope_delay
summed_resp
```

## References

- V. Hohmann (2002), "Frequency analysis and synthesis using a Gammatone filterbank",
  *Acta Acustica united with Acustica* 88(3), 433–442.
- B. R. Glasberg, B. C. J. Moore (1990), "Derivation of auditory filter shapes from
  notched-noise data", *Hearing Research* 47, 103–138.
