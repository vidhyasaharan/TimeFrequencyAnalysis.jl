# Gammatone filterbank analysis (Hohmann 2002) — implementation plan

**Status: planned, not started** (drafted 2026-08-13; revised same day: explicit-bandwidth
construction, long `gammatone_*` names, `filter_coefs` generalised first, per-step
tests + verification plots).

Scope of this first version: the **analysis** side — a bank of 4th-order complex gammatone
filters with ERB-spaced centre frequencies and, by default, ERB-proportional bandwidths (with a
variable scaling factor, and with the bandwidth directly specifiable instead) — plus the
inspection and alignment helpers needed to work with it: impulse responses, frequency
responses, group delay, envelope-peak delay, and delay/phase compensation. Resynthesis (mixer
gains + summation back to a waveform) is **deferred to v2**, but the design keeps that door
open: everything the synthesiser needs (complex subband output, per-channel delays and phase
factors) is already produced by v1.

## 1. What is being built, and why this design

A gammatone filterbank is the standard auditory-motivated time-frequency analysis: each channel
is a bandpass filter whose bandwidth follows the Equivalent Rectangular Bandwidth (ERB) of the
human auditory filter at that centre frequency. Compared with the existing analyses in this
package it fills a distinct niche:

- [`specgram`](../src/spectral_analyses.jl) — uniform frequency resolution, time resolution set
  by the frame length, magnitude only.
- [`periodogram`](../src/spectral_analyses.jl) — arbitrary frequency grid, but still framewise.
- **gammatone filterbank** — ERB frequency resolution, *per-sample* time resolution, and
  **complex** output per channel, separating each subband into envelope (`abs`) and fine
  structure (`angle`) directly.

The chosen design is Hohmann (2002): each channel is a cascade of γ = 4 identical one-pole
**complex** filters. This is preferred over FIR gammatones or the Slaney all-pole real
gammatone because:

1. It is cheap: γ complex multiply–adds per sample per channel, no convolution with long kernels.
2. The complex (analytic-like) output gives the Hilbert envelope for free — no separate
   Hilbert transform pass.
3. It is the basis of an *invertible* analysis–synthesis system (delay + phase compensation +
   mixer + summation), so v2 resynthesis is a natural extension rather than a redesign.
4. Three mature reference implementations exist for cross-validation (below).

### Reference implementations and correspondence

| Concept | Hohmann's MATLAB (Zenodo) | AMT | pyfilterbank | planned TFA |
|---|---|---|---|---|
| single filter design | `gfb_filter_new` | `hohmann2002filter` | `design_filter` | `gammatone_filter` |
| filterbank construction | `gfb_analyzer_new` | `hohmann2002` | `GammatoneFilterbank`, `frequencies_gammatone_bank` | `gammatone_filterbank` |
| run filter/bank over signal | `gfb_filter_process`, `gfb_analyzer_process` | `hohmann2002process` | `fosfilter`, `GammatoneFilterbank.analyze` | `gammatone_filt`, `gammatone_analysis`, `gammatone_cochleagram` |
| frequency response | `gfb_analyzer_zresponse` | `hohmann2002freqz` | `GammatoneFilterbank.freqz` | `filter_resp` / `filter_magresp` methods (via generalised `filter_coefs`) |
| delay + phase compensation | `gfb_delay_new` / `gfb_delay_process` | `hohmann2002delay` | `GammatoneFilterbank.delay` | `gammatone_delay`, `compensate` |
| mixer gains (flatten summed response) | `gfb_mixer_new` | `hohmann2002mixer` | — | *deferred (v2)* |
| resynthesis | `gfb_synthesizer_new` / `_process` | `hohmann2002synth` | — | *deferred (v2)* |

ERB-scale helpers correspond to pyfilterbank's `hertz_to_erbscale` / `erbscale_to_hertz`.

## 2. Mathematics

Notation: sampling rate $f_s$, centre frequency $f_c$, filter order $\gamma$ (default 4),
audiological bandwidth scaling factor $\phi$ (default 1, the user-facing `bw_scale`).

### 2.1 The gammatone filter

The complex-valued analogue gammatone impulse response is

$$g(t) = t^{\gamma-1}\, e^{-2\pi b t}\, e^{i 2\pi f_c t}, \qquad t \ge 0,$$

a complex tone at $f_c$ under a gamma-distribution-shaped envelope. $b$ (in Hz) sets the
envelope decay and therefore the bandwidth. The real gammatone of the classical auditory
literature is the real part; we keep the complex signal because its magnitude is the subband
envelope.

### 2.2 Bandwidth: ERB default, explicit override, scaling factor

Glasberg & Moore's ERB of the auditory filter at $f_c$ (the form used in Hohmann 2002, Eq. 13):

$$\mathrm{ERB}(f_c) = L + \frac{f_c}{Q}, \qquad L = 24.7\ \mathrm{Hz},\quad Q = 9.265.$$

(Equivalently $24.7\,(4.37\, f_c/1000 + 1)$.) The bandwidth actually given to a channel is

$$B(f_c) = \phi \cdot W, \qquad W = \begin{cases}
\mathrm{ERB}(f_c) & \text{default}\\[2pt]
\texttt{bw} & \text{when supplied explicitly,}
\end{cases}$$

i.e. **the bandwidth is an input variable in its own right** (`bw` keyword, per filter; scalar
or per-channel vector for a bank), with the audiological ERB only as its default, and the
scaling factor $\phi$ (`bw_scale`) multiplying whichever is in effect. Bandwidth
(`bw`/`bw_scale`) and channel *density* (filters per ERB, §2.3) remain independent knobs: one
widens each filter, the other packs more of them. A scalar `bw` across a bank gives a
constant-bandwidth complex filterbank — a legitimate non-auditory use of the same machinery.

### 2.3 Centre-frequency grid: the ERB-rate scale

Centre frequencies are spaced uniformly on the ERB-rate (ERB-number) scale

$$E(f) = Q \ln\!\left(1 + \frac{f}{LQ}\right) = 9.265\, \ln\!\left(1 + \frac{f}{228.85}\right),
\qquad f(E) = LQ\left(e^{E/Q} - 1\right),$$

with $LQ = 24.7 \times 9.265 = 228.85$ Hz. Given `fmin`, `fmax`, an anchor `base_frq` and a
density $\nu$ = `filters_per_erb`, the grid is every frequency $f(E_0 + k/\nu)$,
$k \in \mathbb{Z}$, that lies in $[f_{\min}, f_{\max}]$, where $E_0 = E(\texttt{base\_frq})$ —
i.e. steps of $1/\nu$ ERB through a grid point pinned exactly at the anchor. This matches
`gfb_analyzer_new`, whose demo grid (70 Hz – 6.7 kHz, base 1 kHz, $\nu = 1$) gives ~30 channels.

### 2.4 Digital filter design

Each channel is the $\gamma$-fold cascade of one complex one-pole filter:

$$H(z) = k\,\bigl(1 - \tilde a\, z^{-1}\bigr)^{-\gamma}, \qquad \tilde a = \lambda\, e^{i\beta},
\qquad \beta = \frac{2\pi f_c}{f_s},$$

with pole radius $\lambda \in (0,1)$ and normalisation $k$. Implemented as $\gamma$ recursions
of one complex state each ($w_0 \equiv x$):

$$w_j[n] = w_{j-1}[n] + \tilde a\, w_j[n-1], \quad j = 1,\dots,\gamma; \qquad y[n] = k\, w_\gamma[n].$$

Two routes from a bandwidth $B$ to $\lambda$, matching the two call forms of `gfb_filter_new` /
`hohmann2002filter`:

**Route A — ERB matching (the default; used for both the ERB default and an explicit `bw`).**
The analogue prototype's squared magnitude is
$\lvert G(f)\rvert^2 \propto \bigl[1 + ((f-f_c)/b)^2\bigr]^{-\gamma}$, whose own equivalent
rectangular bandwidth is $a_\gamma b$ with (Hohmann 2002, Eq. 14)

$$a_\gamma = \frac{\pi\,(2\gamma-2)!\;2^{-(2\gamma-2)}}{[(\gamma-1)!]^2}
\qquad\Longrightarrow\qquad a_4 = \frac{5\pi}{16} \approx 0.98175.$$

Setting the filter's ERB equal to $B$ and sampling the analogue envelope decay:

$$b = \frac{B(f_c)}{a_\gamma}, \qquad \lambda = e^{-2\pi b / f_s}.$$

**Route B — bandwidth at a given edge attenuation** (gfb-compatible secondary form). Given a
bandwidth $B$ whose edges $f_c \pm B/2$ should sit $A$ dB below the peak (gfb default $A = 3$):
with edge offset $\varphi = \pi B / f_s$ rad/sample and per-stage power ratio
$\alpha = 10^{-A/(10\gamma)}$, requiring
$\lvert H_1 \rvert^2$ at the edge to be $\alpha$ times its peak gives a quadratic
$\lambda^2 + p\lambda + 1 = 0$:

$$p = \frac{-2 + 2\alpha \cos\varphi}{1 - \alpha}, \qquad
\lambda = -\frac{p}{2} - \sqrt{\frac{p^2}{4} - 1},$$

taking the root inside the unit circle (the two roots are reciprocals). This is the
`design_filter` path in pyfilterbank.

The two routes are deliberately different criteria. For $\gamma = 4$ the prototype's −3 dB full
width is $2b\sqrt{2^{1/\gamma}-1} \approx 0.870\,b$, so with Route A the −3 dB width is
$\approx 0.886 \times B$; feeding the same $B$ into Route B pins the −3 dB width to $B$ itself
and yields a filter ≈ 1.13× wider in ERB terms. Route A is the auditory-standard default (as in
gfb/AMT); Route B exists for exact reproduction of gfb's explicit-bandwidth call form.

### 2.5 Normalisation

$$k = 2\,(1-\lambda)^\gamma.$$

Since $\lvert 1 - \tilde a e^{-i\beta}\rvert = 1 - \lambda$, the magnitude response peaks
*exactly* at $f_c$ with gain exactly 2 (the denominator magnitude
$\lvert 1 - \lambda e^{i(\beta - \theta)}\rvert$ is minimised precisely at $\theta = \beta$).
The factor 2 exists because the filter is one-sided in frequency: a real input tone
$A\cos(2\pi f_c t)$ carries amplitude $A/2$ at $+f_c$, so with gain 2 the complex output has
magnitude $\approx A$ — `abs.(y)` is the subband envelope and `real.(y)` the unity-gain
bandpass signal. Both are approximations only insofar as the filter's small residual response
at negative frequencies is neglected; the peak-gain-2 property itself is exact.

### 2.6 Exact discrete impulse response (test anchor, and the `impulse_response` helper)

From the binomial series $(1-x)^{-\gamma} = \sum_n \binom{n+\gamma-1}{\gamma-1} x^n$, the
impulse response of the cascade is **exactly**

$$h[n] = k \binom{n+\gamma-1}{\gamma-1}\, \tilde a^{\,n}, \qquad n \ge 0,$$

for $\gamma = 4$: $h[n] = k\,\frac{(n+1)(n+2)(n+3)}{6}\, \lambda^n e^{i\beta n}$. Two uses:

1. `impulse_response` can be implemented from this closed form **at design time**, before any
   filtering kernel exists — which is what lets impulse responses be plotted as soon as a
   filter can be designed (step 2).
2. It is a machine-precision target for unit tests of the filtering kernel (step 3) — stronger
   than comparing against the sampled analogue $g(t)$, which the digital filter only
   approximates.

Expanding the denominator likewise gives the **polynomial form of the cascade** for the
response machinery (§3.2): $\bigl(1-\tilde a z^{-1}\bigr)^{\gamma} =
\sum_{j=0}^{\gamma} \binom{\gamma}{j} (-\tilde a)^{j} z^{-j}$, i.e. for $\gamma = 4$
denominator coefficients $[1,\ -4\tilde a,\ 6\tilde a^2,\ -4\tilde a^3,\ \tilde a^4]$ and
numerator $[k]$.

### 2.7 Delays: envelope peak vs group delay

Two distinct delay notions, both wanted as helpers:

- **Envelope-peak delay** — the envelope $t^{\gamma-1} e^{-2\pi b t}$ is maximised at
  $$t_{\text{peak}} = \frac{\gamma - 1}{2\pi b}$$
  (equivalently, $\lvert h[n]\rvert \propto \binom{n+3}{3}\lambda^n$ peaks at
  $n \approx -3/\ln\lambda$). This is what the alignment in §2.8 uses.
- **Group delay at $f_c$** — per stage the digital group delay at $\theta = \beta$ is
  $\lambda/(1-\lambda)$ samples, so for the cascade
  $$\tau_g(f_c) = \gamma\,\frac{\lambda}{1-\lambda}\ \text{samples}
  \;\approx\; \frac{\gamma}{2\pi b}\ \text{s},$$
  which is *not* $t_{\text{peak}}$ — the asymmetric envelope puts its maximum earlier than the
  group delay. A `group_delay` helper over an arbitrary frequency grid is computed numerically
  as $-\,d\varphi/d\theta$ from the unwrapped phase response, with the closed form above as its
  value at $f_c$.

Both delays scale as $1/b$: low channels are much slower (see the table in §2.9). This is why
band outputs are misaligned across frequency and compensation is needed before summing
channels or displaying phase-coherent structure.

### 2.8 Delay and phase compensation

Goal: after compensation every channel's response to an impulse peaks — envelope *and* fine
structure — at one common target delay $\tau_d$ (gfb's default: 4 ms), so that
$\sum_k \operatorname{Re}(\tilde y_k[n])$ approximates $\delta[n - n_d]$,
$n_d = \operatorname{round}(\tau_d f_s)$.

Following `gfb_delay_new`: send a unit impulse through the bank and keep only samples
$0 \dots n_d$; in that window find each channel's envelope maximum
$n_{\text{peak}}(k) = \arg\max_{n \le n_d} \lvert h_k[n]\rvert$. Then

- **delay compensation** (integer): $\Delta_k = n_d - n_{\text{peak}}(k) \ \ (\ge 0)$ samples
  of pure delay. Restricting the search window to $[0, n_d]$ makes the clamp implicit: a slow
  low-frequency channel whose true peak lies beyond $n_d$ gets the window edge as its "peak"
  and hence $\Delta_k = 0$ — it cannot be advanced (causality) and remains partially
  misaligned, which the reference design accepts.
- **phase compensation**: multiply the delayed channel by the unit phasor
  $$c_k = \frac{\lvert h_k[n_{\text{peak}}(k)]\rvert}{h_k[n_{\text{peak}}(k)]}, \qquad
  \lvert c_k\rvert = 1,$$
  rotating the channel so it is real and positive at the alignment instant; the channels then
  add constructively at $n_d$.

Compensated subbands: $\tilde y_k[n] = c_k\, y_k[n - \Delta_k]$ (zero-fill at the start, same
length, tail truncated). Note `gfb_delay_new` additionally uses the local phase *slope* at the
envelope peak when computing the phase factor (a refinement related to the peak falling between
samples); the implementation should follow the gfb recipe exactly and validate against
reference outputs rather than re-derive.

### 2.9 Worked numbers (γ = 4, φ = 1, fs = 16 kHz)

| $f_c$ | $\mathrm{ERB}(f_c)$ | $b$ | $\lambda$ | $k$ | $t_{\text{peak}}$ | $\tau_g(f_c)$ |
|---|---|---|---|---|---|---|
| 100 Hz | 35.5 Hz | 36.2 Hz | 0.98590 | 7.9e−8 | 13.2 ms (211 samp) | 17.6 ms |
| 1 kHz | 132.6 Hz | 135.1 Hz | 0.94833 | 1.43e−5 | 3.53 ms (57 samp) | 4.71 ms |
| 4 kHz | 456.4 Hz | 464.9 Hz | 0.83313 | 1.55e−3 | 1.03 ms (16 samp) | 1.37 ms |

Sanity anchors for tests and for the docs page. Note the 100 Hz channel peaks *after* the
default 4 ms alignment target — the clamping case of §2.8 is the norm, not an edge case, at
the bottom of the grid.

## 3. API design

### 3.1 Fit with existing containers and conventions

- The analysis output is a `timefreq` with **complex** components (`S = Complex{T}` is already
  supported): `components` is `nchannels × nsamples`, `frqs` the centre frequencies, `time` the
  per-sample instants `(0:N-1)./fs`, `frames = nothing`. The existing recipes and
  `element_ops` (`amp2db`, …) then work unchanged on the envelope form.
- Every analysis follows the package's dual-output pattern: plain form returns a container,
  `comp()` form returns the bare `Matrix`/`Vector`.
- Precision policy as everywhere else: design arithmetic in `Float64`, but *filtering computes
  in the signal's precision* — a `Float32` signal yields `Complex{Float32}` subbands with the
  scalar coefficients converted at the kernel entry (cheap; they are scalars, not arrays).
- **Frequency responses reuse the existing `filter_coefs` machinery**, generalised to complex
  coefficients first (§3.2, implementation step 1). A `gammatone_filter` converts exactly to a
  `filter_coefs` via the binomial expansion of §2.6, and `filter_resp` / `filter_magresp` /
  `H` / `Hmag` then work as for any other filter. The closed forms (peak gain 2, binomial
  impulse response) remain as independent cross-checks in the tests.
- Filter/response helpers return **plain arrays**, consistent with the existing
  `filter_resp`/`filter_magresp` (which return vectors, not containers). Complex impulse
  responses cannot live in a `signal` container anyway (its samples are real by design).
- No new dependencies: the filtering kernel is a hand-rolled γ-state recursion (`DSP.filt`
  with complex coefficient vectors would work but obscures the per-stage structure);
  `factorial`/`binomial` come from Base (γ is small).

**Memory note**: full-rate complex output is the point of this representation, but it is bulky —
`ComplexF64`, 30 channels, 16 kHz is ~7.7 MB per second of signal (half for `Float32`). v1
accepts this and documents it; a decimated-envelope variant is deferred (§6).

### 3.2 Generalising `filter_coefs` to complex coefficients (step 1)

The current `filter_coefs{T<:AbstractFloat}` holds real `Vector{T}` coefficients, so it cannot
represent the complex-pole gammatone. Generalise following the container pattern already used
by `spectrum`/`timefreq`:

```julia
struct filter_coefs{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    num::Vector{S}
    den::Vector{S}
end
```

with outer constructors deriving `T` (real precision) and `S` (coefficient type) by promotion,
so all existing call sites (`filter_coefs([1],[1,-0.5])`, integer promotion, mixed-precision
promotion) behave identically. `F isa filter_coefs{Float64}` remains true for real-coefficient
filters (`filter_coefs{Float64}` becomes the UnionAll over `S`), so the existing type-check
tests still pass unchanged.

Two required internal fixes discovered by inspection, to be made in the same step:

1. **`H` must not conjugate coefficients.** `H(F, z)` currently evaluates the polynomials with
   `dot(F.num, npz)`, and `dot` conjugates its first argument — harmless for real
   coefficients, wrong for complex ones. Replace with Base's `evalpoly` (ascending-order
   evaluation, exactly the stored layout), which also removes the `powers` allocations.
2. **Normalise the evaluation convention to $z^{-1}$.** `H`/`filter_resp` currently evaluate
   in *positive* powers of `z` (documented quirk: magnitudes standard, phases negated). With
   phase responses and group delay becoming first-class outputs, negated phase would make
   every phase plot and group-delay sign wrong. Switch evaluation to
   `evalpoly(inv(z), coefs)` (on the unit circle `inv(z) == conj(z)`), update the `H`
   docstring, and drop the quirk note. Behaviour change is confined to the *phase* of complex
   responses; magnitudes are identical, and no existing test (or recipe) depends on the phase
   sign — verified by running the full suite.

`Hmag`'s hand-rolled loop generalises by promoting its accumulator from the real coefficient
type to `Complex` of the real precision of `S` and flipping the exponent sign to match the
convention fix.

### 3.3 Frequency-grid additions (`src/frequency_grids.jl`)

```julia
erb(fc)                       #ERB(fc) = 24.7 + fc/9.265, in Hz
freq2erb(f)                   #ERB-rate scale E(f), in ERB numbers  (naming matches freq2θ/freq2z)
erb2freq(E)                   #inverse
erbfreq_array(; fmin=70, fmax=6700, base_frq=1000, filters_per_erb=1)
```

`erbfreq_array` mirrors `logfreq_array`/`linfreq_array`: keyword interface, `Float64` output,
grid pinned at `base_frq` per §2.3. Usable on its own (e.g. as a `periodogram` grid) — worth
exporting regardless of the filterbank. (The short `erb`-family names stay: unlike the analysis
functions they follow the existing `freq2θ`/`logfreq_array` idiom.)

### 3.4 Filter and filterbank types (new file `src/gammatone.jl`)

```julia
struct gammatone_filter{T<:AbstractFloat}
    fs::T; fc::T; order::Int      #γ
    bw::T                         #realised bandwidth B in Hz (after default/override and scaling)
    b::T                          #damping parameter, Hz
    coef::Complex{T}              #ã = λ e^{iβ}
    norm::T                       #k = 2(1-λ)^γ
end

#Primary constructor — Route A. bw overrides the ERB default; bw_scale multiplies either.
gammatone_filter(fs, fc; bw=nothing, bw_scale=1, order=4)
    #B = bw_scale * (bw === nothing ? erb(fc) : bw)

#gfb-compatible secondary form — Route B (edges at fc ± bw/2 sit attenuation_db below peak)
gammatone_filter(fs, fc, bw, attenuation_db; order=4)

struct gammatone_filterbank{T<:AbstractFloat}
    fs::T
    fcs::Vector{T}
    filters::Vector{gammatone_filter{T}}
end
gammatone_filterbank(fs; fmin=70, fmax=min(6700, 0.45fs), base_frq=1000,
                     filters_per_erb=1, order=4, bw=nothing, bw_scale=1)
gammatone_filterbank(fs, fcs; order=4, bw=nothing, bw_scale=1)   #explicit centre frequencies
#fmax caps at 0.45fs so low-fs signals design cleanly with the defaults (decided at step 4)
```

Bank bandwidth semantics: `bw = nothing` → per-channel ERB default; scalar → the same explicit
bandwidth for every channel (constant-bandwidth bank); vector → per-channel bandwidths (length
must match the number of channels). `bw_scale` (scalar) multiplies whichever is in effect.

Storage precision follows `filter_coefs` precedent (promote from the numeric inputs, default
`Float64`); the kernel converts to the signal's precision at call time. Constructor validation:
`0 < fc < fs/2`, `bw > 0`, warn when `fc > 0.45*fs` (the design assumes $f_c \ll f_s/2$; the
pole placement remains valid but the response skirt wraps). Conversion
`filter_coefs(gf::gammatone_filter)` (binomial expansion, §2.6) connects the type to all the
§3.2 response machinery.

### 3.5 Filtering and analyses

```julia
gammatone_filt(gf::gammatone_filter, x::AbstractVector{T})  -> Vector{Complex{T}}
gammatone_filt!(y, gf, x)                                   #in-place, allocation-free
gammatone_filt(fb::gammatone_filterbank, x)                 -> Matrix{Complex{T}}  (nch × N)

gammatone_analysis([comp(),] s::signal, fb)     -> timefreq{T, Complex{T}}   #complex subbands
gammatone_analysis([comp(),] s::signal; <bank kwargs>)                       #designs the bank internally
gammatone_cochleagram([comp(),] s::signal, fb)  -> timefreq{T, T}            #envelope: abs.(subbands) .+ eps(T)
gammatone_cochleagram([comp(),] s::signal; <bank kwargs>)
```

Naming decided: the spelled-out `gammatone_*` forms (`gammatone_analysis` for the complex
subbands, `gammatone_cochleagram` for the envelope representation), paralleling `dft`/`magspec`
in role: complex low-level result vs the magnitude form most users want.
`gammatone_cochleagram` floors with `eps(T)` exactly as `specgram` does, so
`amp2db(gammatone_cochleagram(...))` is always well defined. Titles:
`"Gammatone Filterbank (Hohmann 2002)"` / `"Gammatone Cochleagram"`. Kernel (γ = 4 case;
general γ loops over a small state vector):

```julia
w1 = w2 = w3 = w4 = zero(Complex{T})
@inbounds for n in eachindex(x)
    w1 = muladd(ã, w1, x[n]); w2 = muladd(ã, w2, w1)
    w3 = muladd(ã, w3, w2);   w4 = muladd(ã, w4, w3)
    y[n] = k*w4
end
```

### 3.6 Inspection helpers (the plotting support)

All return plain arrays that plot in one line with the Plots ecosystem; **no new plot recipes
in v1** (cochleagrams already render through the existing `timefreq` recipe).

```julia
impulse_response(gf; dur=nothing)   -> (t, h::Vector{Complex})  #closed form §2.6; dur defaults to decay to −80 dB re peak
filter_coefs(gf)                    -> filter_coefs             #binomial expansion; feeds the machinery below
filter_resp(gf, f)                  -> Vector{Complex}          #= filter_resp(filter_coefs(gf), f, gf.fs)
filter_resp(fb, f)                  -> Matrix                   #nch × nfrq
filter_magresp(gf, f); filter_magresp(fb, f)
summed_resp(fb, f; delay=nothing)   -> Vector{Complex}          #Σ_k responses, optionally compensated
group_delay(gf, f)                  -> Vector                   #numeric −dφ/dθ on the grid, samples
group_delay(gf)                     -> scalar                   #closed form γλ/(1−λ) at fc
envelope_delay(gf)                  -> scalar                   #measured argmax|h|; ≈ (γ−1)/(2πb)
```

Phase responses need no dedicated helper: `angle.(filter_resp(gf, f))` (correct sign after the
§3.2 convention fix; unwrap before differentiating for group delay). The `gf`/`fb` methods of
`filter_resp`/`filter_magresp` are thin wrappers over the generalised machinery — the filter
knows its own `fs`, hence the two-argument signature.

### 3.7 Alignment (delay + phase compensation)

```julia
struct gammatone_delay{T}
    nd::Int                           #target delay in samples
    delays::Vector{Int}               #Δ_k
    phase_factors::Vector{Complex{T}} #c_k
end
gammatone_delay(fb; delay=0.004)                   #measures IRs per §2.8
compensate([comp(),] tf_or_Y, d::gammatone_delay)  -> aligned copy
compensate!(Y, d)
```

Optionally, `gammatone_analysis`/`gammatone_cochleagram` accept
`align::Union{Nothing,Real}=nothing` to fold compensation into the analysis call.

### 3.8 Draft export list

`erb`, `freq2erb`, `erb2freq`, `erbfreq_array`, `gammatone_filter`, `gammatone_filterbank`,
`gammatone_filt`, `gammatone_filt!`, `gammatone_analysis`, `gammatone_cochleagram`,
`impulse_response`, `summed_resp`, `group_delay`, `envelope_delay`, `gammatone_delay`,
`compensate`, `compensate!` (+ new methods on the existing `filter_coefs`, `filter_resp`,
`filter_magresp`). `checkdocs=:exports` in `docs/make.jl` means each needs a docstring
surfaced on a docs page.

## 4. Implementation steps

Each step is a self-contained commit: code + docstrings + tests green **and a visual
verification** before moving on. Test files follow the one-file-per-source-file convention
(`test/gammatone_tests.jl`, plus additions to `filters_tests.jl`, `grids_tests.jl` and
`float32_tests.jl`), registered in `runtests.jl`.

The verification plots accumulate in a dev-only script, **`scripts/gammatone-verification.jl`**
(one clearly labelled section per step, using Plots.jl; not part of the package or test suite).
It keeps every visual check reproducible as the code evolves, and by step 8 it is the ready-made
source of figures and examples for the docs page.

**Step 1 — generalise `filter_coefs` to complex coefficients** (`filters.jl`; §3.2).
Two-parameter struct, promoting constructors, `dot`→`evalpoly`, $z^{-1}$ convention
normalisation, docstring updates.
*Tests*: the **full existing suite passes untouched** (this is the step's headline acceptance
criterion — nothing broken); new tests in `filters_tests.jl`: complex construction and
promotion rules (`isa` checks for `{T}` and `{T,Complex{T}}` forms), a one-pole complex filter
`1/(1 − ã z⁻¹)` has its magnitude peak exactly at the pole angle with gain `1/(1−λ)`, phase of
a pure delay `z⁻¹` at frequency `f` is `−2πf/fs` (convention now standard), real-coefficient
filters keep conjugate-symmetric responses, `Hmag ≈ abs(H)` with complex coefficients.
*Verify visually*: magnitude and phase of the one-pole complex filter across `[0, fs)` — one
asymmetric peak at the pole frequency, no mirror peak (script section 1).

**Step 2 — ERB utilities, filter design, and first response plots**
(`frequency_grids.jl` + `gammatone.jl`: `erb`, `freq2erb`, `erb2freq`, `erbfreq_array`;
`gammatone_filter` with both constructors; `filter_coefs(gf)`; closed-form
`impulse_response`; `filter_resp`/`filter_magresp` methods).
*Tests*: ERB(1 kHz) ≈ 132.63 Hz; `erb2freq(freq2erb(f)) ≈ f` round-trip; grid contains
`base_frq` exactly, spacing uniform in ERB-rate, ~30 channels for the gfb demo parameters;
`0 < λ < 1`; `|filter_resp(gf, [fc])[1]| ≈ 2` exactly (§2.5); Route A: `b ≈ B/a₄` with
`a₄ = 5π/16`; explicit `bw` at fixed `fc` widens/narrows monotonically and
`bw=erb(fc)` reproduces the default; `bw_scale` composes multiplicatively with both; Route B:
magnitude at `fc ± bw/2` is `−attenuation_db` re peak; cross-check coefficients `(ã, k)`
against AMT/pyfilterbank literals (generated offline once, embedded with a provenance
comment); `Float32` inputs produce `gammatone_filter{Float32}`.
*Verify visually* (the step's deliverable the user can look at — script section 2):
```julia
gf = gammatone_filter(16000.0, 1000.0)
t, h = impulse_response(gf)
plot(t, [real.(h) abs.(h)])                       #impulse response + envelope
f = collect(0.0:2.0:8000.0)
plot(f, amp2db.(filter_magresp(gf, f)))           #magnitude response, dB (peak +6 dB at fc)
plot(f, angle.(filter_resp(gf, f)))               #phase response (wrapped)
```
plus overlays: Route A vs Route B at the same `bw`; `bw_scale ∈ {0.5, 1, 2}` family; an
explicit-`bw` filter vs its ERB-default sibling.

**Step 3 — filtering kernel** (`gammatone_filt`, `gammatone_filt!`).
*Tests*: unit impulse in → output matches the exact binomial form of §2.6 to machine precision
(the strongest test in the plan); linearity; steady-state envelope of a tone at `fc` ≈ 1 after
the transient (loose tolerance — negative-frequency leakage); `Float32` signal →
`Complex{Float32}` output with no `Float64` promotion; `gammatone_filt!` allocation-free after
warmup.
*Verify visually*: kernel impulse output overlaid on the closed-form `impulse_response`
(indistinguishable); envelope of a tone burst at `fc` settling to 1; response to the fixture
two-tone signal (script section 3).

**Step 4 — filterbank and the analyses** (`gammatone_filterbank`, bank `gammatone_filt`,
`gammatone_analysis`, `gammatone_cochleagram`, `comp()` forms).
*Tests*: component matrix `nch × N`; `frqs == fb.fcs`; `time == (0:N-1)./fs`; a tone at
`fcs[k]` maximises the mean envelope in channel `k`; scalar-`bw` bank gives every channel the
same `bw`, vector `bw` respected per channel, length mismatch errors;
`gammatone_cochleagram ≈ abs.(gammatone_analysis) .+ eps(T)`; titles; `comp()` forms return
bare arrays; runs on the shared `setup.jl` fixture.
*Verify visually*: all bank magnitude responses overlaid (dB) showing the ERB-proportional
widening — plus the same plot for a constant-`bw` bank; cochleagram heatmap of a synthetic
chirp and of the fixture signal via the existing `timefreq` recipe
(`plot(amp2db(gammatone_cochleagram(sig)))`) next to `specgram` of the same signal (script
section 4).

**Step 5 — delay-related inspection helpers** (`group_delay`, `envelope_delay`,
`summed_resp`).
*Tests*: numeric `group_delay(gf, f)` at `fc` ≈ closed form `γλ/(1−λ)`; `envelope_delay` ≈
`(γ−1)/(2πb)` within one sample; `filter_resp` agrees with the FFT of a long measured impulse
response on the FFT grid; magresp peak location = `fc`.
*Verify visually*: group-delay curves for low/mid/high channels on one axis (the 1/b spread of
§2.9 visible); per-channel envelope-peak times vs the `(γ−1)/(2πb)` curve; summed (
uncompensated) bank response magnitude (script section 5).

**Step 6 — delay and phase compensation** (`gammatone_delay`, `compensate(!)`,
`align` kwarg).
*Tests*: for an impulse input, every channel with `t_peak ≤ τ_d` has its compensated envelope
argmax exactly at `nd`; compensated IR is real-positive at `nd` (phase-factor exactness);
`Σ Re(compensated)` peaks at `nd`; slow channels clamp to `Δ = 0` without error; shape/type
preservation; cross-check `delays`/`phase_factors` against `gfb_delay_new` reference output
for one small bank (embedded literals with provenance).
*Verify visually* (the classic Hohmann figure pair — script section 6): stacked per-channel
`real.(h_k)` impulse responses before vs after compensation (fine structure lining up at 4 ms);
envelope heatmaps before/after; `Σ Re` of the compensated bank ≈ a pulse at `nd`, with
`summed_resp(fb, f; delay=d)` showing the flattened-but-rippled magnitude the v2 mixer will
polish.

**Step 7 — integration and docs.**
`include("gammatone.jl")` after `filters.jl` in the module (it adds methods to names defined
there); exports; module docstring updated; new docs page `docs/src/gammatone.md` (theory
summary = §2 of this doc condensed, usage examples, and figures lifted from
`gammatone-verification.jl`) inserted in `make.jl` pages after "Spectral Analyses"; README
feature bullet; version bump — recommend **0.2.0 → 0.3.0** (a new analysis subsystem plus a
generalised public type is more than the patch-level precedent used for additive plotting
glue).

## 5. Validation strategy

- **Nothing-broken gate at step 1**: the `filter_coefs` generalisation lands only with the full
  pre-existing suite green; the deliberate phase-convention change is confirmed against every
  existing test and documented in the `H` docstring.
- **Closed-form anchors** (preferred wherever available): binomial impulse response (§2.6),
  exact peak gain 2 at `fc` (§2.5), `a₄ = 5π/16`, delay closed forms (§2.7).
- **Cross-implementation literals**: a small set of reference numbers (coefficients, a short IR
  prefix, delays and phase factors for one bank) generated once from AMT (Octave) and/or
  pyfilterbank, embedded in the test files as literals with a comment recording exactly how
  they were produced. No test-time dependency on MATLAB/Python, mirroring how the existing
  tests hardcode known responses.
- **Visual verification per step** via `scripts/gammatone-verification.jl` (§4) — every step
  ends with something to look at, and the script doubles as the docs-figure source.
- **Property tests** on the shared fixture signal: linearity, precision preservation
  (`float32_tests.jl` additions), envelope positivity, channel selectivity.
- **Numerical edge cases**: very low `fc` at high `fs` (λ → 1: verify `Float32` design still
  distinguishes λ from 1 — design in `Float64` then convert), `fc` near the 0.45·fs warning
  threshold, `bw_scale` extremes (0.5, 3), explicit `bw` much narrower/wider than the ERB,
  one-channel banks, signals shorter than `nd`.

## 6. Deferred to v2 (kept in mind, not built now)

- **Mixer + synthesis** (`hohmann2002mixer`/`synth` equivalents): per-channel real gains chosen
  so the delayed-phase-corrected-summed frequency response is ≈ 0 dB at every centre frequency
  (evaluated from the FFT of the summed compensated IRs, iteratively refined), then
  `x̂[n] = Σ_k g_k Re(ỹ_k[n])`. v1's `gammatone_delay` output is exactly the input this needs.
- **Streaming/block processing**: the kernel state is γ complex numbers per channel; a stateful
  processing mode (as in gfb) falls out naturally once there is a use case.
- **Decimated envelope output** (`hop_dur` kwarg producing a framed-rate cochleagram) for
  long-signal memory relief.
- **Real-gammatone convenience** and fractional-delay refinement of the alignment.
- **Makie support** for the new outputs — rides on the existing
  [makie-plot-extension plan](makie-plot-extension.md); nothing gammatone-specific needed.

## 7. Open questions (decide at implementation time)

1. **Route B surface.** Keep the gfb-compatible positional form
   `gammatone_filter(fs, fc, bw, attenuation_db)` as specced, or fold it into the keyword
   constructor (e.g. `attenuation_db=nothing` switching routes)? Current lean: positional, to
   mirror `gfb_filter_new` exactly and keep the primary constructor's `bw` unambiguous
   (Route A).
2. **`gammatone-verification.jl` afterlife**: retire it once the docs page exists, or keep it
   maintained as the manual smoke-test? Current lean: keep — it is cheap and the docs build
   does not exercise plots.

*(Resolved since the first draft: analysis-function naming → long `gammatone_*` forms;
explicit bandwidth input → `bw` keyword on filter and bank, scalar or per-channel, ERB
default, composing with `bw_scale`; `filter_coefs` complex generalisation → yes, done first
as step 1 with machinery reuse.)*

## References

- V. Hohmann (2002), "Frequency analysis and synthesis using a Gammatone filterbank",
  *Acta Acustica united with Acustica* 88(3), 433–442.
- T. Herzke, V. Hohmann (2007), "Improved numerical methods for gammatone filterbank analysis
  and synthesis", *Acta Acustica united with Acustica* 93(3), 498–500. (Numerics used by the
  reference software.)
- Hohmann's original MATLAB (`gfb_*`), citable Zenodo release:
  <https://zenodo.org/records/2643400>
- AMT model pages: <https://www.amtoolbox.org/amt-1.2.0/doc/models/hohmann2002.php> and the
  demo <https://www.amtoolbox.org/amt-1.2.0/doc/demos/demo_hohmann2002.php> (functions
  `hohmann2002`, `hohmann2002filter`, `hohmann2002process`, `hohmann2002delay`,
  `hohmann2002mixer`, `hohmann2002synth`, `hohmann2002freqz`).
- pyfilterbank (faithful Python port incl. delay estimation and phase factors):
  <https://github.com/SiggiGue/pyfilterbank> (docs: <https://siggigue.github.io/pyfilterbank/>).
- B. R. Glasberg, B. C. J. Moore (1990), "Derivation of auditory filter shapes from
  notched-noise data", *Hearing Research* 47, 103–138. (ERB and ERB-rate formulas.)
