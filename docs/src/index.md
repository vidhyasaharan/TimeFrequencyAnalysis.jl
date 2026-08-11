```@meta
CurrentModule = TimeFrequencyAnalysis
```

# TimeFrequencyAnalysis.jl

```@docs
TimeFrequencyAnalysis
```

The package grew out of [SpeechBox.jl](https://github.com/unsw-edu-au/SpeechBox.jl), which
now builds its speech analyses on top of it and re-exports this package's API.

## Design in one paragraph

A signal lives in a [`waveform`](@ref) (samples plus sampling rate). A
[`framed_signal`](@ref) adds a framing geometry without copying anything. Analyses return
their results in [`spectrum`](@ref) (frequency axis) or [`timefreq`](@ref) (time and
frequency axes) containers that keep the provenance of the result — or, when called with
[`comp()`](@ref comp) as the first argument, return the plain component array only. All
containers are parametric in the floating point precision of the signal, and every routine
computes in that precision, so `Float32` in means `Float32` out.

## Quickstart

```julia
using TimeFrequencyAnalysis

fs = 8000.0
x  = cos.(2π*440 .* (1:8000) ./ fs) .+ 0.1 .* white_noise(8000)

s      = waveform(x, fs)               #signal container
frames = framed_signal(s, 0.02, 0.01)  #20 ms frames every 10 ms

spec = magspec(s)                      #spectrum object
tf   = specgram(frames)                #timefreq object (STFT magnitude)
db   = amp2db(tf)                      #elementwise dB, still a timefreq

S    = specgram(comp(), frames)        #the same STFT as a plain Matrix
```

## Contents

```@contents
Pages = ["types.md", "framing.md", "spectral.md", "correlations.md", "filters.md", "plotting.md", "utilities.md", "internals.md"]
Depth = 2
```

## Installation

The package is not yet registered. Install it directly from the repository:

```julia
pkg> add https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl.git
```
