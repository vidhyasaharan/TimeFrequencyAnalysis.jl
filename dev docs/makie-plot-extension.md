# Makie plot extension — design notes (todo, post-migration)

**Status: deferred.** Not part of the SpeechBox → TFA migration (stage 2 of that migration ships the
existing Plots-ecosystem recipes as a RecipesBase package extension). This documents how to add Makie
support when it is picked up, and why it was deferred rather than dropped.

## Why a Makie extension, and why not yet (assessed 2026-08-10)

- **Plots.jl is in maintenance mode**: latest release v1.41.6 (Feb 2026), ~15 commits in the six
  months to Aug 2026, and the v2 rework unreleased after years of development. It still works fine,
  has the largest install base, and is what TFA's existing users use — so the RecipesBase extension
  stays and remains the baseline.
- **Makie is the actively developed ecosystem** (v0.24.13, Jul 2026) and where new users land, but it
  is still 0.x with breaking releases: supporting it implies occasional `[compat]` bumps and
  extension maintenance.
- **Package extensions are the 2026 norm** for library plotting support. DimensionalData.jl ships
  both a RecipesBase and a Makie extension; SignalAnalysis.jl ships a Plots extension. Makie's own
  guidance for library authors: an extension that weak-depends on full `Makie` — never a hard
  dependency, never `MakieCore`, never a backend (GLMakie/CairoMakie).

Deferring keeps the migration surface small; nothing in the migration blocks on, or is blocked by,
this design.

## Design

One package extension, `TimeFrequencyAnalysisMakieExt`, weak-depending on **Makie**. The first
iteration implements `Makie.convert_arguments` methods only — pure data conversion, no custom plot
types — to minimise the surface exposed to Makie's 0.x churn.

### Project.toml changes

```toml
[weakdeps]
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"

[extensions]
TimeFrequencyAnalysisMakieExt = "Makie"

[compat]
Makie = "0.24"
```

Julia 1.10 (TFA's minimum) fully supports package extensions. The RecipesBase extension from
migration stage 2 coexists in the same `[weakdeps]`/`[extensions]` tables.

### ext/TimeFrequencyAnalysisMakieExt.jl — sketch

```julia
module TimeFrequencyAnalysisMakieExt

using TimeFrequencyAnalysis
using Makie

#waveform → line plot of amplitude against time in seconds
Makie.convert_arguments(P::Makie.PointBased, s::waveform) =
    Makie.convert_arguments(P, (0:length(s.x)-1) ./ s.fs, s.x)

#spectrum → line plot of components against frequency in Hz
Makie.convert_arguments(P::Makie.PointBased, sp::spectrum) =
    Makie.convert_arguments(P, sp.frqs, sp.components)

#timefreq → heatmap. components is (nfreq × nframes); Makie expects size (length(x), length(y)),
#so time is x, frqs is y, and the matrix is transposed
Makie.convert_arguments(::Type{<:Makie.Heatmap}, tf::timefreq) =
    (tf.time, tf.frqs, permutedims(tf.components))

end
```

The shapes above are the stable idea; check the exact conversion-trait/type names (`PointBased`,
`Heatmap`, and whether `Makie.expand_dimensions`/plot-type hooks are wanted) against the Makie
version current at implementation time — traits have moved between 0.2x releases.

### Decisions to make at implementation time

1. **Nonuniform frequency grids.** The Plots recipes render `timefreq` in index space with
   `generate_ticks` labels, which displays log-spaced grids (`periodogram` on `logfreq_array`)
   uniformly in log-frequency. Makie options: (a) true coordinates — Makie heatmaps accept
   nonuniform cell-centre vectors natively, so the data is placed correctly but low frequencies
   compress on a linear axis; (b) leave coordinates linear and let users set `yscale = log10` on
   the `Axis`; (c) index space plus manual ticks, mirroring the Plots recipes. Pick one and
   document the choice on the plotting docs page.
2. **Axis labels and titles.** `convert_arguments` converts data only — it cannot set
   "Time (sec)"/"Frequency (Hz)" labels or the container's `title` the way the Plots recipes do.
   Options, cheapest first: leave labelling to the user (fine for a first iteration); export stub
   functions from TFA `src/` (`function plotspecgram end`, methods added in the extension — the
   pattern Makie's docs recommend) that build a fully labelled `Axis`; or a full Makie `@recipe`
   if recipe-level axis control has landed in the Makie version in use.
3. **No silent rescaling.** Conversions must not apply dB scaling implicitly; users write
   `heatmap(amp2db(tf))` explicitly, matching the Plots-recipe behaviour.

### Testing

Makie loads headless without a backend, and `convert_arguments` is unit-testable without rendering:
add `Makie` to `test/Project.toml` and assert the returned tuples (time axis values, the
transposition, eltype preservation for `Float32` signals). No image-comparison tests needed. Note
that Makie is a heavy test dependency — either accept the CI cost consciously or gate these tests
in their own test group.

### Versioning

Purely additive → patch bump under 0.x semver (SpeechBox's `[compat] TimeFrequencyAnalysis`
tolerates patch releases). Update SpeechBox's compat only if this lands together with a minor bump.

## References

- Makie library-author guidance: <https://docs.makie.org/stable/explanations/recipes>
- Worked weakdep-extension example: <https://github.com/jkrumbiegel/MakiePkgExtTest>
- Dual-ecosystem precedent (RecipesBase + Makie extensions): <https://github.com/rafaqz/DimensionalData.jl>
- Plots-extension precedent in DSP: <https://github.com/org-arl/SignalAnalysis.jl>
