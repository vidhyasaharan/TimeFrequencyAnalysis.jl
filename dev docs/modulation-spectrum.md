# Modulation spectrum analysis (design + build plan)

**Status: planned 2026-08-17. Not started.** Target release: **v0.4.0**. Branch: **`main`**.

Adds the generic half of an amplitude modulation spectrum analysis: from a complex gammatone
subband representation to *power as a function of carrier frequency × modulation rate*. The
sonar-specific half (octave bands, beams, band-limiting) lives downstream in
[SONAR.jl](https://github.com/vidhyasaharan/SONAR.jl) — see `SONAR/dev docs/phase2-modulation-spectrum.md`.

Driven by a SONAR requirement, but nothing here knows about beams or octaves: this is the same
representation used for speech (Atlas & Shamma's joint acoustic/modulation frequency),
psychoacoustics (Ewert & Dau's envelope power spectrum model, Houtgast & Steeneken's MTF) and
machine condition monitoring (envelope spectra).

## Working method

Build **one function at a time**, each complete only when it has an implementation, a docstring
in the house style, tests, and inspection before the next starts. [Build order](#build-order) is
a dependency order.

**Test both precisions.** `test/float32_tests.jl` gains a modulation block; individual test
files spot-check `Float32` inline as they already do.

## Scope

| In | Out (later) |
| --- | --- |
| `modfreq` container (carrier × modulation rate) | Modulation *filterbank* (Ewert & Dau style bandpass bank on the envelope) |
| `power_envelope` — subband envelope, decimated | Time-varying modulation spectra (3-D: carrier × rate × time) |
| `modulation_spectrum` — Welch over the envelope | Cyclostationary / spectral-coherence estimators |
| `welch_psd` — a general PSD estimator, the missing primitive | Synthesis / inverse |
| `resample` for matrices, and a `Float32` fix for the existing method | Modulation-domain filtering |
| Delay-compensation helpers: `envelope_delay(fb)`, `default_align`, `compensation_lead`, `compensate(...; trim)` | |
| Plot recipe for `modfreq` | |
| `scripts/modulation-verification.jl` — visual checks per build step | |

## The chain

```
signal
  ↓ per-channel complex gammatone filt!        (exists)
  ↓ optional delay/phase compensation + drop the zero-filled lead
  ↓ abs2 (or abs)                              → envelope
  ↓ anti-alias + integer decimation            → envelope @ fs_env
  ↓ Welch: frame → dc → window → rfft → |·|²   → averaged over frames
                                               → modfreq
```

## Design decisions

### `modfreq` does not carry the originating signal

`spectrum` and `timefreq` both anchor to the `signal` they came from, because their axes index
*into* it. A modulation spectrum is a statistic *over* a whole span — nothing in it points back
at a sample. Carrying the signal would also pin ~90 MB per beam-hour downstream, and SONAR's
`fullband_modspec` has **four** originating signals from four different physical arrays, so the
field could only ever be `nothing` in the most-used case.

It carries `bws` (per-channel analysis bandwidth) instead of the filterbank, because that is the
one thing needed downstream, it concatenates when a bank does not, and it states the per-row
validity limit: column `j` of row `k` is only meaningful while `mod_frqs[j] ≲ bws[k]/2`.

### Never materialise the full complex subband matrix

`filt(fb, x)` allocates `n_channels × n_samples` `Complex{T}` — 96 MB for a 60 s SONAR span at
32 carriers, **5.8 GB for an hour**. `power_envelope` loops channels with the in-place
single-filter kernel `filt!(ybuf, fb.filters[k], x)`, reusing one complex buffer, and decimates
each row before storing it. Peak memory becomes one channel's buffer plus the decimated output.

### `max_mod_frq` is primary, `fs_env` is the override

`max_mod_frq` is the physically meaningful quantity — a user knows they care about modulation to
125 Hz, not about a sample rate. `M = max(1, floor(Int, fs/(2.5·max_mod_frq)))`, `fs_env = fs/M`.
The 2.5× (not 2×) keeps `max_mod_frq` clear of the anti-alias transition band. Integer decimation
only: `rate = 1//M` is exact, single-phase and phase-linear, and `fs_env` need not be round.

Anti-aliasing uses DSP's own Kaiser-window resample filter with taps **converted to `T`**. This
matters: `DSP.resample_filter` returns `Float64` taps and the output promotes, so the existing
`resample(s::signal, fs_new)` silently returns `Float64` for a `Float32` signal today. That is a
latent bug in v0.3.0 and is fixed as part of this work.

### Delay compensation is on by default, and the zero-fill is dropped

`compensate!` shifts row `k` right by `Δ_k` and **zero-fills its first `Δ_k` samples**. Those
fabricated zeros must not reach an analysis, so when compensation is applied the first
`maximum(d.delays)` columns are dropped from every row (before decimation), keeping the matrix
rectangular and the channels time-aligned. **The time axis keeps its original origin** — it
starts at `compensation_lead(d)/fs`, not 0 — so the result still refers to the same instants as
the input.

The existing fixed 4 ms default target is wrong for low-frequency work: a channel at
`fc = 77 Hz` with `B = 66 Hz` has an envelope delay of ~7.4 ms, so it never peaks inside a 4 ms
search window and is left partially misaligned (a failure mode `gammatone_delay`'s docstring
documents and accepts). `default_align(fb) = (maximum(envelope_delay(fb)) + 1)/fb.fs` derives the
target from the bank itself, so no channel is left behind. `align = :auto` uses it.

Per-channel `|DFT|²` is shift-invariant, so compensation does not move a peak *within* one row's
modulation spectrum. It matters for anything comparing channels, which is why it is on by default
and why the helpers are general-purpose and live in `gammatone.jl`, not here.

### DC handling: four options, not two

Subtracting the envelope mean and dividing by it are **independent** operations. In the DFT,
subtraction changes only bin 1 while division rescales every bin, so the four combinations are
genuinely distinct and all are exposed:

| `dc` | Operation | Units | DC bin |
|---|---|---|---|
| `:keep` | `x` | physical | mean carrier power |
| `:remove` | `x − μ` | physical | ~0 |
| `:divide` | `x / μ` | dimensionless | 1 |
| `:index` (default) | `(x − μ)/μ` | dimensionless | ~0 |

`:index` is the true modulation index and the default. Published DEMON removes DC; the
MTF/modulation-depth tradition normalises by the mean; neither alone is sufficient. Implemented
as two internal flags resolved from the symbol.

### `welch_psd` — what gets reused, and what cannot

TFA has no PSD estimator today. That is a real gap, useful well beyond modulation work, so this
is a new exported function rather than a private helper — but it should reuse as much of the
existing machinery as it can.

**Reused as-is:**

- the **`spectrum{T,T}` container** — a Welch estimate is exactly a vector of components on a
  frequency grid, which is what `spectrum` already is. No new type.
- **`framed_signal` + `view_frame` + `number_signal_frames`** for framing. `view_frame` returns
  views rather than copies, and `num_signal_frames` already excludes the zero-padded tail frames
  that would bias an average downward.
- the **`window`** function, one `plan_rfft`, and `rfftfreq` — i.e. `specgram`'s whole skeleton
  (`spectral_analyses.jl:86-108`).

**Cannot be reused:**

- **`enframe`.** Two reasons, neither of them the shift bug any more — that was found while
  designing this and is **fixed** in `00b01e8`, which reroutes `enframe` through
  `framed_signal`/`view_frame` so the requested shift is now honoured exactly. What remains:

  1. It returns **all** `num_frames` columns, including the trailing ones that extend past the
     end of the signal and are zero padded. Averaging a zero-padded frame into a power estimate
     biases it downward. Welch must use `num_signal_frames` only, which is what `view_frame`
     lets it do directly.
  2. It materialises the whole `frame_size × num_frames` matrix. In `modulation_spectrum` that
     is one such allocation *per carrier channel* — 116 of them at SONAR's defaults — on top of
     the envelope it already holds.

  Since `enframe` is now itself a thin wrapper over `framed_signal` + `view_frame`, using those
  two directly is not bypassing it so much as calling the same machinery one level down, without
  the materialisation.
- **`specgram`.** Returns *magnitude*, floors with `eps(T)`, has no per-frame DC handling and no
  power scaling, and keeps every frame when only their mean is wanted.
- **`magspec` / `dft`.** Single transform, no averaging, no power scaling; `dft` allocates a
  fresh buffer and window on every call.

So `welch_psd` shares `specgram`'s structure but accumulates `|X|²` into a running mean.

### The estimator

For frames ``x_f[n]``, `f = 1…F`, of length `N`, window `w`, `ndft = M`, and `x̃` the frame after
the `dc` rule:

```
X_f[k] = Σ_{n=0}^{N-1} w[n] x̃_f[n] e^{-i2πkn/M}

P[k]   = (c · α_k / F) Σ_{f=1}^{F} |X_f[k]|²

α_k    = 1 for k = 0 and k = M/2, else 2          (one-sided fold)
c      = 1/(Σ_n w[n])²       for scaling = :spectrum
c      = 1/(f_s Σ_n w[n]²)   for scaling = :density
```

`:spectrum` (default) reads a sinusoid of amplitude `A` as `A²/2` at its bin — right for the
discrete shaft/blade lines this exists to find. `:density` gives `Σ_k P[k]·Δf ≈ var(x)` — right
for noise floors, and frame-length independent. Both properties are asserted in the tests.

Why average at all: a single periodogram of noise has ~100 % standard error **however long the
record**. Only averaging `F` frames reduces it, by ``1/\sqrt{F}``. That is the whole reason this
function exists rather than `magspec`.

### `welch_psd` vs `specgram`, for the docstring

| | `specgram` | `welch_psd` |
|---|---|---|
| Returns | `timefreq` — one column per frame | `spectrum` — one vector, frames averaged |
| Quantity | magnitude `\|X\|` | power `\|X\|²`, scaled |
| Frames | all kept (time axis) | averaged away (variance ∝ 1/F) |
| Floor | `+ eps(T)` on every bin | none |
| DC / scaling | neither | `dc` and `scaling` keywords |
| Answers | *when* something happened | *how much* power at each frequency |

### `envelope = :power | :amplitude`

The literature is genuinely split: power in the MTF/EPSM and classical square-law DEMON
tradition, magnitude in recent underwater-ML gammatone work. The switch is one operator, so
expose it. Default `:power`.

## Defaults

At `fs = 3125` (the SONAR case), the useful settings:

| `max_mod_frq` | M | `fs_env` | mod Nyquist | `frame_dur` | ndft | Δf_mod |
|---|---|---|---|---|---|---|
| **125 (default)** | 10 | 312.5 Hz | 156.25 Hz | 10.0 s | **3125** (=5⁵) | **0.10 Hz** |
| 50 | 25 | 125 Hz | 62.5 Hz | 4.0 s | 500 | 0.25 Hz |
| 10 | 125 | 25 Hz | 12.5 Hz | 10.0 s | 250 | 0.10 Hz |

The 125 Hz default is set by the physics: blade rate reaches ~140 Hz for a 20 Hz shaft on a
7-bladed propeller, and the widest channel of a SONAR bank (`B = 315 Hz` at `fc = 1232 Hz`) can
only carry modulation to ~158 Hz anyway. `frame_dur = 10.0` matches published sonar integration
times (5–80 s).

## The `B ≳ 2·f_m` rule

A filter of bandwidth `B` can only resolve a modulation at `f_m` if both sidebands `fc ± f_m`
sit inside it. This is the constraint that binds at the bottom of a low-frequency bank, and it
is what `bws` exists to expose.

**No citation was found for it** in a literature search — it is correct physics but
bibliographically unsupported, so it is documented as *derived, not cited*, and pinned by a real
test (a `B/f_m ∈ {0.5, 1, 2, 4}` sweep) rather than a comment. Atlas & Shamma (2003) and Ewert &
Dau (2000, JASA 108:1181–1196) are the likely sources; neither was accessible.

## Build order

Steps 3–6 each carry a **"Verify visually"** item. Following the precedent of
`scripts/gammatone-verification.jl`, these live in a companion script
**`scripts/modulation-verification.jl`**, one section per step, run with
`julia> include("scripts/modulation-verification.jl")` (Plots is not a package dependency). It is
not part of the test suite; by step 8 it doubles as the source of figures for the docs page.

Tests assert numbers; these plots catch the things a passing assertion hides — an axis
transposed, a decade of scaling error, an envelope that tracks but lags, a peak in the right bin
sitting on a floor that should not be there.

1. **Delay helpers** — `envelope_delay(fb)`, `default_align(fb)`, `gammatone_delay(fb; delay = :auto)`,
   `compensate(...; trim = false)`, `compensation_lead(d)`, in `src/gammatone.jl`. Tests in
   `test/gammatone_tests.jl`. *Independent of everything below; can land on its own.* **DONE.**
2. **`resample` for matrices** + the `Float32` tap fix on the existing method, in
   `src/sampling.jl` via a shared private `_resample_taps`. Tests in `sampling_tests.jl` and
   `float32_tests.jl`.
3. **`welch_psd`**, in `src/spectral_analyses.jl`. Tests: Parseval under `:density`, sinusoid
   amplitude under `:spectrum`, and the `dc` orthogonality identities.

   **Verify visually:**
   - *Variance reduction.* One 60 s white-noise record, PSD estimated with `F = 1, 4, 16, 64`
     frames, overlaid on a log-power axis with the true `σ²/f_Nyq` as a horizontal line. The
     scatter must visibly shrink as ``1/\sqrt{F}`` while the mean stays put. This is the single
     most informative plot for the whole function — it shows both that the scaling is right (the
     curves sit on the line, not a decade off) and that averaging is doing what it is for.
   - *Scaling conventions.* `A cos(2πf₀t) + noise` for `A = 1`, plotted twice: under `:spectrum`
     the line reads `A²/2 = 0.5` at its bin regardless of frame length; under `:density` it grows
     with frame length while the *noise floor* stays fixed. Annotate both with the expected value.
     Getting these two backwards is the classic PSD bug and this plot makes it obvious.
   - *The four `dc` modes.* One modulated envelope, four PSDs overlaid. `:keep`/`:remove` must be
     indistinguishable everywhere except bin 1; likewise `:divide`/`:index`; and the two
     dimensionless modes must be a constant factor below the physical pair. Plot the
     bin-1-excluded ratio to make the orthogonality visible rather than asserted.
   - *Window and overlap.* The effective noise bandwidth and the frame-to-frame overlap actually
     used, printed alongside the requested values — the guard against the `enframe` shift trap
     described above, applied to whatever framing `welch_psd` ends up using.
4. **`modfreq`** type in `src/types.jl` + `element_ops.jl` methods. Tests in `types_tests.jl`,
   `element_ops_tests.jl`.

   **Verify visually:**
   - *Axis orientation, with a deliberately asymmetric pattern.* Fill `components` with a known
     ramp that differs along the two axes (e.g. `10*i + j`), give it a log carrier grid and a
     linear modulation grid, and plot. Carrier must be on **y**, modulation rate on **x**, and
     the value at a named `(fc, f_m)` must land where the axes say. A transposed heatmap is the
     easiest mistake to make here and the hardest to notice later, because a plausible-looking
     image survives it.
   - *Tick labels on a log carrier axis.* Ticks must read as frequencies (78, 110, 156 …), not
     row indices, and stay legible with 32 rows — the recipe works in index space, so this is
     where `generate_ticks` gets checked.
   - *`pow2db` round-trip.* The same object plotted linear and in dB, confirming the element-op
     methods return a `modfreq` with axes intact rather than a bare matrix.
5. **`power_envelope`**, new `src/modulation.jl`. Tests: decimation length and unit DC gain,
   anti-alias rejection, the gain-2 amplitude identity, the bandwidth-visibility sweep, and the
   compensation/zero-lead behaviour.

   **Verify visually:**
   - *Envelope tracks the waveform.* For a single AM tone, plot `real(y_k)` and `abs(y_k)` on one
     time axis for one channel, over a couple of modulation periods. The envelope must ride the
     peaks — not lag them, not sit at half height. Add the true `A(1 + m cos 2πf_m t)` as a dashed
     line; any gain error shows immediately.
   - *Decimation is honest.* The envelope spectrum before and after decimation, overlaid up to
     the new Nyquist. They must coincide in the passband and the decimated one must roll off into
     the anti-alias transition — no folded-back energy. Repeat with a deliberately out-of-band
     modulation (100 Hz at `fs_env = 125`) and mark where its alias *would* appear at 25 Hz; the
     plot must show nothing there.
   - *Alignment.* A broadband click through the bank, plotted as an image (channels × time) over
     a short window, with `align = nothing` and `align = :auto` side by side. Uncompensated shows
     the characteristic diagonal sweep of channel delays; compensated shows one vertical ridge.
     Mark `compensation_lead` with a vertical line and confirm the trimmed output starts after it.
   - *The `B ≳ 2f_m` rule.* An image of modulation depth recovered versus `B/f_m` swept over
     `{0.5, 1, 2, 4, 8}`. This is the plot that stands in for the citation the literature search
     could not find, so it belongs in the docs page too.
6. **`modulation_spectrum`**, same file. Tests: AM tone lands in the right cell, peak level
   `2m²A⁴`, second harmonic 18 dB down at `2f_m`, two carriers at two rates with no cross-talk,
   all four `dc` modes.

   **Verify visually:**
   - *The headline figure — two carriers, two rates.* 150 Hz modulated at 3 Hz plus 600 Hz
     modulated at 11 Hz, as a `modfreq` heatmap in dB. Two isolated blobs at the two `(carrier,
     rate)` cells, and nothing at the two off-diagonal cells. This one image demonstrates the
     entire point of a two-dimensional representation and should open the docs page.
   - *Harmonic structure.* A single carrier at `m = 0.5`, one carrier row plotted as a line
     spectrum, with `f_m` and `2f_m` marked and the predicted 18 dB drop annotated. Sweep `m` over
     `{0.1, 0.25, 0.5, 1.0}` on the same axes to show the harmonic growing as `m` does.
   - *Dynamic range and floor.* The same signal with noise added at a few SNRs, to see where the
     line disappears into the floor — the practical question a user actually has.
   - *Against a spectrogram, on the same data.* `specgram` beside `modulation_spectrum` for a
     signal with AM sidebands: the sideband spacing visible in the spectrogram must equal the
     modulation rate the modulation spectrum reports. This cross-check is what will be run on real
     SONAR beams at acceptance, so it is worth having working on synthetic data first.
7. **Plot recipe** in `ext/`. Tests in `recipes_tests.jl`. *(The step-4 and step-6 plots above
   assume it; until it lands they run against `heatmap(ms.components)` directly.)*
8. **Finish**: `float32_tests.jl` block; `docs/src/modulation.md` + `make.jl`; `modfreq` in
   `types.md`; `welch_psd` in `spectral.md`; `default_align`/`compensation_lead` in `gammatone.md`;
   `modfreq` bullet in `plotting.md`; module docstring; version 0.4.0.

`makedocs` runs `checkdocs = :exports`, so **every new export must appear in a docs page or the
build fails**.

## Acceptance

- `Pkg.test()` green, `docs/make.jl` green.
- A 200 Hz carrier modulated at 7 Hz peaks at exactly that `(carrier, rate)` cell, at the level
  theory predicts.
- `Float32` in → `modfreq{Float32,Float32}` out, matching `Float64` to `rtol = 1e-3`.
- A broadband click peaks at the same output index in every channel with `align = :auto`, and at
  different indices with `align = nothing`.
- `scripts/modulation-verification.jl` runs end to end and every figure reads as its build-order
  entry describes — in particular the two-carrier/two-rate heatmap shows two isolated blobs and
  nothing at the off-diagonal cells, and the `B/f_m` sweep shows modulation depth appearing
  around `B/f_m ≈ 2`.

## Related work

This chain is a **multiband DEMON** (Detection of Envelope Modulation On Noise). Published
precedent for the multi-carrier form: *Multiband Enhancement for DEMON Processing Algorithms*,
and patent EP2009459B1 ("improved DEMON analysis using sub-band signals"). The rigorous
formulation of the same carrier × modulation map is the **cyclostationary** one (spectral/cyclic
coherence) — see Marmara, *Time-cyclic frequency representation of propeller tonals*,
[10.1007/s11760-024-03575-6](https://doi.org/10.1007/s11760-024-03575-6). This filter-bank
version is an approximation of it, with the carrier axis chosen by filterbank rather than by
integration band; a cyclic-coherence estimator would be a natural v0.5 addition.
