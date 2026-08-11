# TimeFrequencyAnalysis

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/dev/)
[![Build Status](https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)

Time-frequency analysis of uniformly sampled signals in Julia:

- **Containers** — `signal` (samples + sampling rate), `framed_signal` (framing geometry),
  `spectrum` and `timefreq` (analysis results carrying their frequency/time axes)
- **Framing** — `enframe`, `view_frame`, `extract_frame`, `frame_energy`
- **Spectral analyses** — `dft`, `magspec`, `specgram` (STFT magnitude), `periodogram` on
  arbitrary (e.g. log-spaced) frequency grids
- **Utilities** — window functions, frequency grids, resampling, synthetic signal generators,
  elementwise `log`/`log10`/`amp2db`/`pow2db` on analysis results

All types are parametric in the floating point precision of the signal (`Float64`,
`Float32`, ...) and every routine computes in that precision.

The package originated as the time-frequency core of
[SpeechBox.jl](https://github.com/vidhyasaharan/SpeechBox.jl), which now builds on it and
re-exports its API.

## Installation

Not yet registered; install from the repository:

```julia
pkg> add https://github.com/vidhyasaharan/TimeFrequencyAnalysis.jl.git
```

## Quickstart

```julia
using TimeFrequencyAnalysis

fs = 8000.0
x  = cos.(2π*440 .* (1:8000) ./ fs) .+ 0.1 .* white_noise(8000)

s      = signal(x, fs)
frames = framed_signal(s, 0.02, 0.01)   # 20 ms frames every 10 ms
tf     = amp2db(specgram(frames))       # STFT magnitude spectrogram in dB
```

See the [documentation](https://vidhyasaharan.github.io/TimeFrequencyAnalysis.jl/dev/) for
the full API.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
