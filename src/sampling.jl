#Helpers relating durations in seconds to sample counts, and resampling.
#
#Resampling wraps the polyphase resampler from DSP.jl, adding two things DSP does not do on
#its own: the filter taps are built at the precision of the data (see _resample_taps), and an
#exactly rational rate is used whenever the two sampling rates permit one (see _resample_rate).

"""
    time2nsamples(dur, fs)

Number of samples spanned by a duration of `dur` seconds at sampling rate `fs` Hz, rounded
to the nearest integer.
"""
time2nsamples(dur::Real, fs::Real) = Int(round(dur*fs))


## Resampling

#Number of polyphase branches used when the rate cannot be written as a small fraction; this
#is DSP's own default for the arbitrary-rate resampler
const _RESAMPLE_PHASES = 32

#Largest numerator or denominator accepted for an exactly rational rate. Beyond this the
#rational resampler would need more polyphase branches than the arbitrary-rate one, so the
#arbitrary path is the cheaper approximation
const _MAX_RESAMPLE_RATIO = 1000

#Express fs_new/fs as an exact fraction when a small one exists, and as a Float64 otherwise.
#Which one comes back decides the algorithm: a Rational rate resamples with one polyphase
#branch per output sample, exactly and phase linearly, whereas a Float rate goes through the
#arbitrary-rate resampler, which interpolates between _RESAMPLE_PHASES branches. Integer
#decimation - the case the modulation envelope work needs - always lands on the exact path.
function _resample_rate(fs::Real, fs_new::Real)
    (fs > 0) || error("The current sampling rate must be positive, got $fs")
    (fs_new > 0) || error("The new sampling rate must be positive, got $fs_new")
    r = convert(Float64,fs_new)/convert(Float64,fs)
    q = rationalize(r; tol = 1e-9)
    if(max(numerator(q),denominator(q)) ≤ _MAX_RESAMPLE_RATIO)
        return q
    end
    return r
end

#Resampling taps at the precision of the data. DSP always designs in Float64 and the polyphase
#output promotes to match its taps, so resampling a Float32 signal with the default taps would
#return Float64 - the one place in the package where precision did not propagate.
_resample_taps(::Type{T}, rate::Union{Integer,Rational}, rel_bw::Real, attenuation::Real) where {T<:AbstractFloat} =
    convert(Vector{T},DSP.Filters.resample_filter(rate, rel_bw, attenuation))

_resample_taps(::Type{T}, rate::AbstractFloat, rel_bw::Real, attenuation::Real) where {T<:AbstractFloat} =
    convert(Vector{T},DSP.Filters.resample_filter(rate, _RESAMPLE_PHASES, rel_bw, attenuation))

#DSP takes the phase count as a trailing positional argument on the arbitrary-rate methods only
_apply_resample(x::AbstractVector, rate::Union{Integer,Rational}, h::Vector) = DSP.Filters.resample(x, rate, h)
_apply_resample(x::AbstractVector, rate::AbstractFloat, h::Vector) = DSP.Filters.resample(x, rate, h, _RESAMPLE_PHASES)
_apply_resample(x::AbstractArray, rate::Union{Integer,Rational}, h::Vector, dims::Int) = DSP.Filters.resample(x, rate, h; dims)
_apply_resample(x::AbstractArray, rate::AbstractFloat, h::Vector, dims::Int) = DSP.Filters.resample(x, rate, h, _RESAMPLE_PHASES; dims)

"""
    resample(s::signal, fs_new [; <keyword arguments>])
    resample(x::AbstractVector, fs, fs_new [; <keyword arguments>])
    resample(X::AbstractArray, fs, fs_new [; dims = 2, <keyword arguments>])

Resample to the new sampling rate `fs_new`, returning a new [`signal`](@ref) for a `signal`
input and a plain array otherwise. For an array of more than one dimension, `dims` selects
the axis to resample along; the default `dims = 2` resamples each **row**, which is the
layout of the channels-by-samples matrices produced by
[`gammatone_analysis`](@ref) and [`power_envelope`](@ref).

Wraps the polyphase resampler from
[`DSP.jl`](https://docs.juliadsp.org/stable/contents/), which lowpass filters before
decimating (or after interpolating), so the result is free of aliasing without any further
filtering by the caller. The anti-alias filter is a Kaiser-windowed FIR of unit DC gain, so
the mean of the resampled data equals the mean of the input.

The output length is `ceil(length(x)*fs_new/fs)`.

**Precision follows the data.** The taps are designed in `Float64` and converted to the
element type of the input, so a `Float32` input resamples in `Float32` and returns `Float32`.

**The rate decides the algorithm.** When `fs_new/fs` is an exact fraction with numerator and
denominator no larger than $(_MAX_RESAMPLE_RATIO) — which covers integer decimation, integer
interpolation, and every simple ratio between them — the exact rational resampler is used, one
polyphase branch per output sample. Otherwise the arbitrary-rate resampler interpolates between
$(_RESAMPLE_PHASES) branches, and `fs_new/fs` is met only approximately.

### Keyword Arguments
- `dims` : Axis to resample along, for the array forms [Default = 2, i.e. along rows]
- `rel_bw` : Relative bandwidth of the anti-alias filter [Default = 1]
- `attenuation` : Stopband attenuation of the anti-alias filter in dB [Default = 60]

Resampling to the rate the data already has is a copy: there is nothing to filter, and the
resampler cannot design a filter for a rate of 1 in any case (its cutoff would sit exactly at
Nyquist).

# Throws
- `ErrorException` if either sampling rate is not positive.
"""
function resample(s::signal{T}, fs_new::Number; rel_bw::Real = 1, attenuation::Real = 60) where {T<:AbstractFloat}
    rate = _resample_rate(s.fs, fs_new)
    (rate == 1) && return signal(copy(s.x), fs_new)
    h = _resample_taps(T, rate, rel_bw, attenuation)
    return signal(_apply_resample(s.x, rate, h), fs_new)
end

function resample(x::AbstractVector{T}, fs::Real, fs_new::Real; rel_bw::Real = 1, attenuation::Real = 60) where {T<:AbstractFloat}
    rate = _resample_rate(fs, fs_new)
    (rate == 1) && return copy(x)
    h = _resample_taps(T, rate, rel_bw, attenuation)
    return _apply_resample(x, rate, h)
end

function resample(X::AbstractArray{T}, fs::Real, fs_new::Real; dims::Int = 2, rel_bw::Real = 1, attenuation::Real = 60) where {T<:AbstractFloat}
    (1 ≤ dims ≤ ndims(X)) || error("dims = $dims is not an axis of a $(ndims(X))-dimensional array")
    rate = _resample_rate(fs, fs_new)
    (rate == 1) && return copy(X)
    h = _resample_taps(T, rate, rel_bw, attenuation)
    return _apply_resample(X, rate, h, dims)
end
