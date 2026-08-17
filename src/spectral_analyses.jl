#Spectral and spectro-temporal analyses: DFT magnitude spectra, STFT spectrograms and
#periodograms evaluated on arbitrary frequency grids.
#
#Each analysis has two forms selected by dispatch (see the comp marker in types.jl):
#the plain form returns a spectrum/timefreq container; the comp() form returns only the
#component array and is the one to use in inner loops.

"""
    dft(x [, <keyword argument>])

Compute the complex DFT coefficients of the array `x` (positive frequencies only, computed
with a real-input FFT). The signal is windowed before transforming.

### Keyword Arguments
- `ndft` : Number of DFT points [Default = `length(x)`]. If `ndft` is less than the length of `x`, the number of DFT points is raised to `length(x)`
- `wtype` : Window type to use [Default = "hanning"], see [`window`](@ref)
"""
function dft(x::AbstractVector{T}; ndft::Int = length(x), wtype::String="hanning") where {T<:AbstractFloat}
    len = length(x)
    win = window(T,len;wtype)
    buflen = max(ndft,len)
    buf = zeros(T, buflen)
    buf[1:len] = x.*win
    return rfft(buf)
end


"""
    magspec(s [, <keyword arguments>])
    magspec(x, fs = 2π [, <keyword arguments>])

Compute the DFT magnitude spectrum of the [`signal`](@ref) `s` (or of the signal in
array `x` with sampling rate `fs`). Returns a [`spectrum`](@ref) object.

    magspec(comp(), x, fs [, <keyword arguments>])

As above, but return the magnitude components and the frequency grid as a tuple of plain
arrays.

    magspec(comp(), x [, <keyword arguments>])

Return only the magnitude components of the DFT of `x`.

### Keyword Arguments
- `ndft` : Number of DFT points [Default is the length of the signal array]
- `wtype` : Window type to use [Default is "hanning"], see [`window`](@ref)
"""
magspec(s::signal;ndft::Int = length(s.x), wtype::String="hanning") = magspec(s.x, s.fs; ndft, wtype)

function magspec(x::AbstractVector{<:AbstractFloat},fs::Real=2π; ndft::Int = length(x), wtype::String="hanning")
    mspec,frqs = magspec(comp(), x, fs; ndft, wtype)
    return spectrum(signal(x,fs), mspec, frqs, "DFT Magnitude Spectrum")
end

function magspec(::comp, x::AbstractVector{<:AbstractFloat},fs::Real; ndft::Int = length(x), wtype::String="hanning")
    ndft = max(ndft, length(x)) #dft() raises the number of DFT points to the signal length; keep the frequency grid consistent with that
    mspec = magspec(comp(), x; ndft, wtype)
    frqs = rfftfreq(ndft, fs)
    return mspec,frqs
end

magspec(::comp, x::AbstractVector{<:AbstractFloat}; ndft::Int = length(x), wtype::String = "hanning") = abs.(dft(x; ndft, wtype))


"""
    specgram(frames::framed_signal [; <keyword arguments>])
    specgram(comp(), frames::framed_signal [; <keyword arguments>])
    specgram([comp(),] x, fs [; <keyword arguments>])

Compute the short-time DFT magnitude spectrogram of the [`framed_signal`](@ref) `frames`
(or of the signal in array `x` with sampling rate `fs`, framed with the durations given by
the keyword arguments). Only the `num_signal_frames` frames that lie entirely within the
signal are analysed. With the [`comp()`](@ref comp) argument the output is the plain matrix
of magnitude components (one column per frame); otherwise it is a [`timefreq`](@ref) object
whose frequency axis matches the rows of the component matrix.

A floor of `eps(T)` is added to every component so that logarithms of the spectrogram are
always well defined.

### Keyword Arguments
- `wtype` : Window type [Default is "hanning"], see [`window`](@ref)
- `ndft` : Number of DFT points [Default is the optimal FFT length that is at least the frame length]
- `frame_dur` : Frame duration in secs, for the signal-array forms [Default = 0.02]
- `frame_shift_dur` : Interval between consecutive frames in secs, for the signal-array forms [Default = 0.01]
"""
function specgram(::comp, sig_frames::framed_signal{T}; ndft::Int = nextfastfft(sig_frames.frame_length), wtype::String = "hanning") where {T<:AbstractFloat}
    len = sig_frames.frame_length
    nfft = max(len,ndft) #Never use fewer DFT points than the frame length
    nframes = sig_frames.num_signal_frames

    win = window(T,len;wtype=wtype)

    buf = zeros(T,nfft) #Zero padded buffer holding one windowed frame at a time
    rfp = plan_rfft(buf) #Precomputed real-input FFT plan, reused for every frame
    nrfft = length(rfp*buf) #Number of (positive) frequency components
    mspec = Matrix{T}(undef,nrfft,nframes)
    for i ∈ 1:nframes
        frame = view_frame(sig_frames,i)
        buf[1:1:len] = win.*frame #Window, then store in the padded buffer
        mspec[:,i] = abs.(rfp*buf)
    end

    #Floor the spectrogram so a log spectrogram is always well defined
    for i ∈ eachindex(mspec)
        mspec[i] += eps(T)
    end
    return mspec
end

function specgram(sig_frames::framed_signal; ndft::Int = nextfastfft(sig_frames.frame_length), wtype::String="hanning")
    mspec = specgram(comp(), sig_frames; ndft, wtype)
    nfft = max(sig_frames.frame_length, ndft) #Number of DFT points used by the computation above
    frqs = rfftfreq(nfft,sig_frames.signal.fs)
    return timefreq(sig_frames,mspec,frqs)
end

function specgram(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; frame_dur::Real=0.02, frame_shift_dur::Real=0.01, ndft::Int = nextfastfft(time2nsamples(frame_dur,fs)), wtype::String="hanning")
    sig_frames = framed_signal(x,fs,frame_dur,frame_shift_dur)
    return specgram(comp(), sig_frames;ndft, wtype)
end

function specgram(x::Array{<:AbstractFloat},fs::Number;frame_dur::Real=0.02, frame_shift_dur::Real=0.01, ndft::Int = nextfastfft(time2nsamples(frame_dur,fs)), wtype::String="hanning")
    sig_frames = framed_signal(x,fs,frame_dur,frame_shift_dur)
    return specgram(sig_frames; ndft, wtype)
end


## Welch power spectrum
#
#specgram keeps every frame and returns magnitude; welch_psd averages the frames away and
#returns scaled power. The two share their skeleton (one window, one rfft plan, view_frame over
#num_signal_frames) but not their purpose - see the welch_psd docstring.
#
#Framing goes through framed_signal/view_frame, the same machinery specgram uses. enframe is not
#used even though it now honours its frame shift exactly: it returns all num_frames columns
#including the zero padded trailing ones, which would bias a power average downward, and it
#materialises the whole frame matrix, which the per-channel loop in modulation_spectrum would
#pay for once per carrier.

#A dc symbol resolves to two independent flags. Subtracting the frame mean changes only the DC
#bin; dividing by it rescales every bin. All four combinations are meaningful, so both are kept
#rather than bundling them into one "normalise" switch.
function _dc_flags(dc::Symbol)
    (dc === :keep)   && return (false, false)
    (dc === :remove) && return (true,  false)
    (dc === :divide) && return (false, true)
    (dc === :index)  && return (true,  true)
    error("dc must be one of :keep, :remove, :divide or :index, got :$dc")
end

#Apply the dc rule to one frame, window it, and write it into the zero padded transform buffer.
#The buffer beyond the frame length is never written, so it stays zero from allocation.
function _prepare_frame!(buf::AbstractVector{T}, frame::AbstractVector{T}, win::AbstractVector{T},
                         subtract::Bool, normalise::Bool) where {T<:AbstractFloat}
    len = length(frame)
    if(!(subtract | normalise))
        @inbounds for n ∈ 1:len
            buf[n] = win[n]*frame[n]
        end
        return buf
    end

    μ = sum(frame)/len
    if(normalise)
        #Dividing by the mean is only meaningful when the frame has one to speak of. Compare
        #against the frame's own RMS rather than against zero: a mean far below the signal scale
        #would otherwise pass and produce an enormous, meaningless modulation index.
        rms = sqrt(sum(abs2,frame)/len)
        if(rms == 0)
            @inbounds for n ∈ 1:len
                buf[n] = zero(T)
            end
            return buf                  #an all-zero frame contributes no power under any rule
        end
        (abs(μ) > sqrt(eps(T))*rms) || error("dc = :divide and dc = :index scale each frame by its mean, but a frame here has a mean of $μ against an RMS of $rms; use dc = :remove or dc = :keep for a signal that is not offset from zero")
        iμ = one(T)/T(μ)
        if(subtract)
            @inbounds for n ∈ 1:len
                buf[n] = win[n]*(frame[n]*iμ - one(T))
            end
        else
            @inbounds for n ∈ 1:len
                buf[n] = win[n]*frame[n]*iμ
            end
        end
    else
        m = T(μ)
        @inbounds for n ∈ 1:len
            buf[n] = win[n]*(frame[n] - m)
        end
    end
    return buf
end

#Sum |X|^2 over the whole frames of a framed_signal into P, reusing one buffer, one transform
#and one output vector. Only num_signal_frames is walked, so view_frame always hands back a view
#and never one of its zero padded tail frames.
#
#modulation_spectrum runs this once per carrier channel. Since signal wraps a Vector it cannot
#hold a row of an envelope matrix directly, so that caller keeps ONE row buffer and ONE
#framed_signal around it and refills the buffer per channel: no allocation per channel, and the
#frame reads stay contiguous rather than striding across the matrix.
function _welch_accumulate!(P::AbstractVector{T}, buf::AbstractVector{T}, X::AbstractVector{Complex{T}},
                            rfp, win::AbstractVector{T}, frames::framed_signal{T},
                            subtract::Bool, normalise::Bool) where {T<:AbstractFloat}
    fill!(P, zero(T))
    fbuf = view(buf,1:frames.frame_length)
    for i ∈ 1:frames.num_signal_frames
        _prepare_frame!(fbuf, view_frame(frames,i), win, subtract, normalise)
        mul!(X, rfp, buf)
        @inbounds for k ∈ eachindex(P)
            P[k] += abs2(X[k])
        end
    end
    return P
end

#Average over frames, apply the scaling convention, and fold the negative frequencies onto the
#positive ones. DC has no negative counterpart, and neither does Nyquist when it exists (only
#for an even transform length), so those two bins are not doubled.
function _welch_finish!(P::AbstractVector{T}, win::AbstractVector{T}, fs::Real, nfft::Int,
                        nframes::Int, scaling::Symbol) where {T<:AbstractFloat}
    c = if(scaling === :spectrum)
            one(T)/T(sum(win))^2
        else
            one(T)/(T(fs)*T(sum(abs2,win)))
        end
    c /= T(nframes)
    has_nyquist = iseven(nfft)
    n = length(P)
    @inbounds for k ∈ 1:n
        α = (k == 1 || (has_nyquist && k == n)) ? one(T) : T(2)
        P[k] *= α*c
    end
    return P
end

"""
    welch_psd([comp(),] frames::framed_signal [; <keyword arguments>])
    welch_psd([comp(),] s::signal [; <keyword arguments>])
    welch_psd([comp(),] x, fs [; <keyword arguments>])

Welch estimate of the power spectrum of a signal: the signal is cut into overlapping frames,
each frame is windowed and transformed, and the resulting power spectra are **averaged**.
Returns a [`spectrum`](@ref); with [`comp()`](@ref comp) the plain component vector for the
`framed_signal` form, or a `(components, frqs)` tuple for the `signal` and array forms (the
same split as [`specgram`](@ref) and [`magspec`](@ref)). Only the `num_signal_frames` frames
that lie entirely within the signal are used, so no zero padded frame biases the average.

### The estimate

For frames ``x_f[n]``, ``f = 1…F``, of length `N`, window ``w``, and `ndft = M`, with
``\\tilde{x}_f`` the frame after the `dc` rule below:

```math
X_f[k] = \\sum_{n=0}^{N-1} w[n]\\, \\tilde{x}_f[n]\\, e^{-i2\\pi kn/M}
\\qquad
P[k] = \\frac{c\\,\\alpha_k}{F} \\sum_{f=1}^{F} |X_f[k]|^2
```

``\\alpha_k = 1`` at ``k = 0`` and at ``k = M/2`` (which exists only for even `M`) and
``\\alpha_k = 2`` elsewhere, folding each negative frequency onto its positive counterpart, and
``c`` is set by `scaling`:

- `scaling = :spectrum` (default): ``c = 1/(\\sum_n w[n])^2``. A sinusoid of amplitude ``A``
  reads ``A^2/2`` at its bin, **independent of frame length** — the convention to use when
  looking for discrete lines.
- `scaling = :density`: ``c = 1/(f_s \\sum_n w[n]^2)``, a power spectral density in units per
  Hz, satisfying ``\\sum_k P[k]\\,\\Delta f ≈ \\mathrm{var}(x)``. Independent of frame length
  for noise — the convention to use when characterising a noise floor.

Averaging is the whole point. A single periodogram of noise has a standard error of about 100%
of its own value **however long the record**; only averaging ``F`` frames reduces it, by
``1/\\sqrt{F}``. Lengthening the frames buys frequency resolution, not reliability.

### Removing the mean

`dc` selects what happens to each frame before it is windowed. Subtracting the mean and
dividing by it are independent operations, so all four combinations are available:

| `dc` | frame becomes | units | DC bin |
|:--- |:--- |:--- |:--- |
| `:keep` | ``x`` | physical | mean power |
| `:remove` (default) | ``x - \\mu`` | physical | ~0 |
| `:divide` | ``x/\\mu`` | dimensionless | 1 |
| `:index` | ``(x-\\mu)/\\mu`` | dimensionless | ~0 |

The two operations behave quite differently in the result. **Dividing** rescales every bin by
exactly ``1/\\mu^2``, whatever the window: `:divide` is `:keep` divided by ``\\mu^2``, and
`:index` is `:remove` divided by ``\\mu^2``, to machine precision. **Subtracting** removes the
frame's constant component, which reaches only the DC bin under `wtype = "rect"` but spreads
over the mainlobe and sidelobes of any tapered window — with a Hann window and a mean-2.5
envelope, the bin next to DC carries 3.15 under `:keep` against 6.4e-10 under `:remove`.
Suppressing that leakage is the reason to subtract the mean before windowing rather than
trusting the window to contain it, and it is why `:remove` rather than `:keep` is the default.

The mean is taken per frame, not over the whole signal, so a slow drift in level is removed
along with the constant offset.

`:divide` and `:index` are meant for strictly positive data such as the envelopes behind
[`modulation_spectrum`](@ref), where `:index` is the modulation index and is that function's
default. They error on a frame whose mean is negligible against its own RMS, and an all-zero
frame contributes no power rather than erroring. That check is a numerical guard against
dividing by nothing, not a judgement about whether the data suits the mode: a signal that
merely happens to straddle zero will divide by a small mean and return large, unstable numbers
without complaint. For a general signal — which may well be centred on zero — `:remove` is the
right choice and is the default here.

### Difference from [`specgram`](@ref)

| | `specgram` | `welch_psd` |
|:--- |:--- |:--- |
| Returns | [`timefreq`](@ref), one column per frame | [`spectrum`](@ref), frames averaged |
| Quantity | magnitude ``\\|X\\|`` | power ``\\|X\\|^2``, scaled |
| Frames | all kept, forming a time axis | averaged away, variance ∝ ``1/F`` |
| Floor | `eps(T)` added to every bin | none |
| `dc`, `scaling` | not offered | both |
| Answers | *when* something happened | *how much* power at each frequency |

### Keyword Arguments
- `frame_dur` : Frame duration in secs, for the `signal` and array forms [Default = 0.02]
- `frame_shift_dur` : Interval between consecutive frames in secs [Default = `frame_dur/2`, i.e. 50% overlap]
- `ndft` : Number of DFT points [Default is the optimal FFT length that is at least the frame length]
- `wtype` : Window type [Default = "hanning"], see [`window`](@ref)
- `dc` : Frame mean handling, one of `:keep`, `:remove`, `:divide`, `:index` [Default = `:remove`]
- `scaling` : `:spectrum` or `:density` [Default = `:spectrum`]

# Throws
- `ErrorException` if the signal is shorter than one frame, if `dc` or `scaling` is not one of
  the values above, or if `:divide`/`:index` meets a frame whose mean is negligible against its
  own RMS.
"""
function welch_psd(::comp, frames::framed_signal{T}; ndft::Int = nextfastfft(frames.frame_length),
                   wtype::String = "hanning", dc::Symbol = :remove, scaling::Symbol = :spectrum) where {T<:AbstractFloat}
    len = frames.frame_length
    nfft = max(len,ndft) #Never use fewer DFT points than the frame length
    nframes = frames.num_signal_frames
    (nframes ≥ 1) || error("The signal holds $(length(frames.signal.x)) samples, fewer than the $len of one frame, so there is nothing to average")
    (scaling ∈ (:spectrum,:density)) || error("scaling must be :spectrum or :density, got :$scaling")
    subtract,normalise = _dc_flags(dc)

    win = window(T,len;wtype)
    buf = zeros(T,nfft)  #Zero padded buffer holding one windowed frame at a time
    rfp = plan_rfft(buf) #Precomputed real-input FFT plan, reused for every frame
    nrfft = div(nfft,2) + 1
    X = Vector{Complex{T}}(undef,nrfft)
    P = Vector{T}(undef,nrfft)

    _welch_accumulate!(P,buf,X,rfp,win,frames,subtract,normalise)
    return _welch_finish!(P,win,frames.signal.fs,nfft,nframes,scaling)
end

function welch_psd(frames::framed_signal; ndft::Int = nextfastfft(frames.frame_length),
                   wtype::String = "hanning", dc::Symbol = :remove, scaling::Symbol = :spectrum)
    P = welch_psd(comp(), frames; ndft, wtype, dc, scaling)
    nfft = max(frames.frame_length,ndft)
    frqs = rfftfreq(nfft,frames.signal.fs)
    title = (scaling === :density) ? "Welch Power Spectral Density" : "Welch Power Spectrum"
    return spectrum(frames.signal,P,frqs,title)
end

function welch_psd(s::signal; frame_dur::Real = 0.02, frame_shift_dur::Real = frame_dur/2,
                   ndft::Int = nextfastfft(time2nsamples(frame_dur,s.fs)), wtype::String = "hanning",
                   dc::Symbol = :remove, scaling::Symbol = :spectrum)
    return welch_psd(framed_signal(s,frame_dur,frame_shift_dur); ndft, wtype, dc, scaling)
end

function welch_psd(::comp, s::signal; frame_dur::Real = 0.02, frame_shift_dur::Real = frame_dur/2,
                   ndft::Int = nextfastfft(time2nsamples(frame_dur,s.fs)), wtype::String = "hanning",
                   dc::Symbol = :remove, scaling::Symbol = :spectrum)
    frames = framed_signal(s,frame_dur,frame_shift_dur)
    P = welch_psd(comp(), frames; ndft, wtype, dc, scaling)
    return P, rfftfreq(max(frames.frame_length,ndft),s.fs)
end

welch_psd(x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = welch_psd(signal(x,fs); kwargs...)

welch_psd(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = welch_psd(comp(), signal(x,fs); kwargs...)


## Complex exponential helpers used by the periodogram

"""
    cexp(T, f, fs, N)
    cexp(f, fs, N)

Unit-norm complex negative exponential of length `N` at frequency `f` Hz for sampling rate
`fs`, with element type `Complex{T}` (default `Complex{Float64}`). The periodogram is the
squared magnitude of inner products with these sequences.
"""
function cexp(::Type{T}, f::Real, fs::Real, N::Int) where {T<:AbstractFloat}
    k = T(2π)*(T(f)/T(fs))
    return (1/sqrt(T(N)))*exp.((-im*k).*(1:N))
end

cexp(f::Real,fs::Real,N::Int) = cexp(Float64,f,fs,N)

"""
    cexp_proj_matrix(T, frqs, fs, N)
    cexp_proj_matrix(frqs, fs, N)

Matrix whose row `i` is [`cexp`](@ref)`(T, frqs[i], fs, N)`; multiplying a signal by it
projects the signal onto all the complex exponentials at once.
"""
function cexp_proj_matrix(::Type{T}, frqs::AbstractVector{<:Real},fs::Real,N::Int) where {T<:AbstractFloat}
    nfrqs = length(frqs)
    proj_matrix = zeros(Complex{T},nfrqs,N)
    for i in eachindex(frqs)
        proj_matrix[i,:] = cexp(T,frqs[i],fs,N)
    end
    return proj_matrix
end

cexp_proj_matrix(frqs::AbstractVector{<:Real},fs::Real,N::Int) = cexp_proj_matrix(Float64,frqs,fs,N)


"""
    periodogram([comp(),] x, fs, frqs[; wtype="hanning"])
    periodogram([comp(),] x, fs [; wtype="hanning"[, fmin=10[, fmax=fs/2]]])
    periodogram([comp(),] sig_frames::framed_signal, frqs[; wtype="hanning"])
    periodogram([comp(),] sig_frames::framed_signal[; wtype="hanning"[, fmin=10[, fmax=sig_frames.signal.fs/2]]])

Compute the periodogram of the signal in array `x` with sampling rate `fs` at the
frequencies in `frqs` — or at frequencies equally spaced on the log scale between `fmin` and
`fmax` (see [`logfreq_array`](@ref)) — as the squared magnitude of the inner product between
the windowed signal and a unit-norm complex exponential at each frequency. Unlike the
FFT-based [`magspec`](@ref)/[`specgram`](@ref), the frequency grid is arbitrary.

When the input is a [`framed_signal`](@ref), the periodogram of each frame is computed.
With the [`comp()`](@ref comp) argument the components are returned as a plain
`Vector`/`Matrix`; otherwise a [`spectrum`](@ref)/[`timefreq`](@ref) object is returned.
"""
function periodogram(x::AbstractVector{<:AbstractFloat},fs::Number,frqs::AbstractVector{<:Number};wtype::String="hanning")
    components = periodogram(comp(), x,fs,frqs;wtype)
    s = signal(x,fs)
    return spectrum(s,components,frqs,"Periodogram")
end

function periodogram(x::AbstractVector{<:AbstractFloat},fs::Number;wtype::String="hanning",fmin::Number=10,fmax::Number=fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(x,fs,frqs;wtype=wtype)
end

function periodogram(sig_frames::framed_signal,frqs ;wtype::String="hanning")
    pspec = periodogram(comp(), sig_frames, frqs; wtype)
    return timefreq(sig_frames,pspec,frqs)
end

function periodogram(sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(sig_frames, frqs, wtype = wtype)
end

function periodogram(::comp, x::AbstractVector{T}, fs::Real, frqs::AbstractVector{<:Real}; wtype::String="hanning") where {T<:AbstractFloat}
    len = length(x)
    win = window(T,len;wtype=wtype)
    ip = x.*win
    nfrqs = length(frqs)
    proj = zeros(T,nfrqs)
    for i ∈ eachindex(proj)
        proj[i] = abs2(dot(ip,cexp(T,frqs[i],fs,len)))
    end
    return proj
end

function periodogram(::comp, x::AbstractVector{<:AbstractFloat},fs::Real;wtype::String="hanning",fmin::Real=10,fmax::Real=fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(comp(), x,fs,frqs;wtype=wtype)
end

function periodogram(::comp, sig_frames::framed_signal{T}, frqs ;wtype::String="hanning") where {T<:AbstractFloat}
    len = sig_frames.frame_length
    nframes = sig_frames.num_signal_frames
    fs = sig_frames.signal.fs
    pspec = Matrix{T}(undef,length(frqs),nframes)
    nfrqs = length(frqs)
    win = window(T,len;wtype=wtype)

    #The full projection matrix is precomputed once and reused for every frame; the
    #transpose is collected so that the inner loop walks down contiguous columns
    ce_array = collect(transpose(cexp_proj_matrix(T,frqs,fs,len)))

    for i ∈ 1:nframes
        ip = extract_frame(sig_frames,i)
        ip .*= win
        @views for j in 1:nfrqs
            pspec[j,i] = abs2(dot(ip,ce_array[:,j]))
        end
    end
    return pspec
end

function periodogram(::comp, sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(comp(), sig_frames, frqs, wtype = wtype)
end
