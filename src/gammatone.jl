#Gammatone filters after Hohmann (2002): each filter is a cascade of `order` identical
#one-pole complex filters H(z) = k/(1 - ã z⁻¹)^order with ã = λe^{iβ}, giving a complex
#(analytic-like) subband output whose magnitude is the band envelope. This file implements
#the single-filter design and its inspection helpers; the filterbank, cochleagram analyses
#and delay/phase alignment build on it (see "dev docs/gammatone-filterbank.md"). Frequency
#responses are evaluated through the filter_coefs machinery (filters.jl) via the exact
#polynomial expansion of the cascade.

"""
    gammatone_filter{T<:AbstractFloat}
    gammatone_filter(fs, fc; bw = nothing, bw_scale = 1, order = 4)
    gammatone_filter(fs, fc, bw, attenuation_db; order = 4)

A complex gammatone filter of order `order` (default 4) at centre frequency `fc` Hz for
sampling rate `fs` Hz, after Hohmann (2002): a cascade of `order` one-pole complex filters
with combined transfer function ``H(z) = k(1 - ãz^{-1})^{-order}`` and pole
``ã = λe^{iβ}``, ``β = 2πf_c/f_s``. The magnitude response peaks at `fc` with gain exactly
2 (``k = 2(1-λ)^{order}``), so that for a real input the complex output carries the
amplitude of the input band: its magnitude is the band envelope and its real part the
bandpass signal.

The keyword constructor sets the filter's equivalent rectangular bandwidth to
`bw_scale*bw` Hz, where `bw` defaults to the auditory bandwidth [`erb`](@ref)`(fc)`, via
the bandwidth-to-damping relation of Hohmann (2002) Eqs. 13–17 (see
[`a_gamma`](@ref TimeFrequencyAnalysis.a_gamma)). The 4-argument form instead places the
band edges `fc ± bw/2` exactly `attenuation_db` dB below the peak (Hohmann 2002,
Eqs. 11–12 — the explicit-bandwidth form of the reference implementations, whose default
attenuation is 3 dB). The two conventions differ: for a 4th-order filter the −3 dB width
is ≈ 0.886 × the equivalent rectangular bandwidth.

# Fields
- `fs::T` : sampling rate in Hz
- `fc::T` : centre frequency in Hz
- `order::Int` : filter order (number of cascaded one-pole sections)
- `bw::T` : design bandwidth in Hz (after default/override and scaling)
- `b::T` : envelope damping in Hz (the impulse-response envelope decays as ``e^{-2πbt}``)
- `coef::Complex{T}` : the pole ``ã = λe^{iβ}``
- `norm::T` : the normalisation ``k = 2(1-λ)^{order}``
"""
struct gammatone_filter{T<:AbstractFloat}
    fs::T
    fc::T
    order::Int
    bw::T
    b::T
    coef::Complex{T}
    norm::T
end

function gammatone_filter(fs::Real, fc::Real; bw::Union{Nothing,Real} = nothing, bw_scale::Real = 1, order::Int = 4)
    T = float(promote_type(typeof(fs),typeof(fc)))
    fs64,fc64 = check_gammatone_args(fs,fc,order)
    W = (bw === nothing) ? erb(fc64) : convert(Float64,bw)
    B = convert(Float64,bw_scale)*W
    (B > 0) || error("Filter bandwidth must be positive")
    b = B/a_gamma(order)
    λ = exp(-2π*b/fs64)
    return build_gammatone_filter(T,fs64,fc64,order,B,b,λ)
end

function gammatone_filter(fs::Real, fc::Real, bw::Real, attenuation_db::Real; order::Int = 4)
    T = float(promote_type(typeof(fs),typeof(fc)))
    fs64,fc64 = check_gammatone_args(fs,fc,order)
    B = convert(Float64,bw)
    A = convert(Float64,attenuation_db)
    (B > 0) || error("Filter bandwidth must be positive")
    (A > 0) || error("attenuation_db must be positive (dB below the peak)")
    φ = π*B/fs64
    α = 10.0^(-A/(10*order))
    p = (-2 + 2*α*cos(φ))/(1 - α)
    (p^2/4 - 1 ≥ 0) || error("The bandwidth/attenuation combination is not realisable at this sampling rate")
    λ = -p/2 - sqrt(p^2/4 - 1) #root of λ² + pλ + 1 inside the unit circle
    b = -fs64*log(λ)/(2π) #damping implied by the pole radius, so the delay relations hold for both routes
    return build_gammatone_filter(T,fs64,fc64,order,B,b,λ)
end

"""
    check_gammatone_args(fs, fc, order)

Validate gammatone design parameters and return `(fs, fc)` converted to `Float64`, the
design arithmetic precision. Errors for out-of-range values; warns when `fc > 0.45fs`,
where the one-sided design's neglect of the response wrap-around starts to bite.
"""
function check_gammatone_args(fs::Real, fc::Real, order::Int)
    fs64 = convert(Float64,fs)
    fc64 = convert(Float64,fc)
    (fs64 > 0) || error("Sampling rate must be positive")
    (0 < fc64 < fs64/2) || error("Centre frequency must lie in (0, fs/2)")
    (1 ≤ order ≤ 10) || error("Filter order must lie in 1:10")
    (fc64 > 0.45*fs64) && @warn "Centre frequency is close to the Nyquist frequency; the gammatone design assumes fc ≪ fs/2 and the response skirt will wrap"
    return fs64,fc64
end

"""
    a_gamma(order)

Bandwidth factor ``a_γ`` of Hohmann (2002) Eq. 14, relating the damping `b` of an
`order`-th order gammatone filter to its equivalent rectangular bandwidth ``a_γ b``. For
the default order 4, ``a_4 = 5π/16 ≈ 0.98175``.
"""
a_gamma(order::Int) = π*factorial(2*order-2)*2.0^(-(2*order-2))/factorial(order-1)^2

#Shared constructor tail: from the design values (Float64) to a filter stored at precision T
function build_gammatone_filter(::Type{T}, fs::Float64, fc::Float64, order::Int, B::Float64, b::Float64, λ::Float64) where {T<:AbstractFloat}
    β = 2π*fc/fs
    ã = λ*cis(β)
    k = 2*(1-λ)^order
    return gammatone_filter{T}(convert(T,fs),convert(T,fc),order,convert(T,B),convert(T,b),convert(Complex{T},ã),convert(T,k))
end


"""
    filter_coefs(gf::gammatone_filter)

Exact rational form of the gammatone cascade: numerator `[k]` and denominator the binomial
expansion of ``(1 - ãz^{-1})^{order}``, as a complex-coefficient [`filter_coefs`](@ref).
Connects the gammatone filter to the response machinery ([`filter_resp`](@ref),
[`filter_magresp`](@ref)).
"""
function filter_coefs(gf::gammatone_filter{T}) where {T<:AbstractFloat}
    γ = gf.order
    den = Vector{Complex{T}}(undef,γ+1)
    for j ∈ 0:γ
        den[j+1] = binomial(γ,j)*(-gf.coef)^j
    end
    return filter_coefs{T,Complex{T}}([Complex{T}(gf.norm)], den)
end

"""
    filter_resp(gf::gammatone_filter, f)

Complex frequency response of the gammatone filter at each analogue frequency of the
vector `f` in Hz (the sampling rate comes from the filter itself), evaluated through
[`filter_coefs`](@ref)`(gf)`.
"""
filter_resp(gf::gammatone_filter, f::AbstractVector{<:Real}) = filter_resp(filter_coefs(gf),f,gf.fs)

"""
    filter_magresp(gf::gammatone_filter, f)

Magnitude response of the gammatone filter at each analogue frequency of the vector `f` in
Hz (the sampling rate comes from the filter itself).
"""
filter_magresp(gf::gammatone_filter, f::AbstractVector{<:Real}) = filter_magresp(filter_coefs(gf),f,gf.fs)


"""
    impulse_response(gf::gammatone_filter; dur = nothing)

Impulse response of the gammatone filter as a tuple `(t, h)`: sample times `t` in seconds
and the complex response `h`, computed from the exact closed form
``h[n] = k\\binom{n+order-1}{order-1}ã^{\\,n}`` (no filtering involved, so the impulse
response is available directly from the design). `abs.(h)` is the response envelope and
`real.(h)` the real gammatone impulse response. The duration defaults to the time the
envelope takes to decay to −80 dB relative to its peak; pass `dur` in seconds to override.
"""
function impulse_response(gf::gammatone_filter{T}; dur::Union{Nothing,Real} = nothing) where {T<:AbstractFloat}
    N = (dur === nothing) ? default_ir_length(gf) : max(1,Int(round(convert(Float64,dur)*convert(Float64,gf.fs))))
    γ = gf.order
    #Computed in Float64 (the design precision) and stored at the filter's precision: the
    #recursion is a long product chain that would shed digits in Float32
    ã = convert(Complex{Float64},gf.coef)
    k = convert(Float64,gf.norm)
    h = Vector{Complex{T}}(undef,N)
    c = 1.0                    #binomial(n+γ-1, γ-1), updated recursively
    an = one(Complex{Float64}) #ãⁿ
    h[1] = k
    for n ∈ 1:N-1
        c *= (n+γ-1)/n
        an *= ã
        h[n+1] = convert(Complex{T},k*c*an)
    end
    t = convert(Vector{T},(0:N-1)./convert(Float64,gf.fs))
    return t,h
end

"""
    default_ir_length(gf)

Number of samples after which the impulse-response envelope of `gf` has decayed to −80 dB
relative to its peak: the default length used by [`impulse_response`](@ref). Found by
walking the envelope ratio ``|h[n+1]|/|h[n]| = λ(n+order)/(n+1)`` outwards from the
envelope peak.
"""
function default_ir_length(gf::gammatone_filter)
    γ = gf.order
    λ = convert(Float64,abs(gf.coef))
    n = max(1,Int(round(-(γ-1)/log(λ)))) #envelope peak position
    r = 1.0
    while (r > 1e-4) && (n < 10_000_000) #-80 dB in amplitude, with a runaway guard
        r *= λ*(n+γ)/(n+1)
        n += 1
    end
    return n+1
end


## Filtering kernel

"""
    gammatone_filt(gf::gammatone_filter, x)

Filter the real signal array `x` with the gammatone filter `gf` and return the complex
subband signal: `abs.` of the output is the band envelope and `real.` the
bandpass-filtered signal (see [`gammatone_filter`](@ref) for why the gain-2 normalisation
makes both carry the input band's amplitude). The cascade is applied as its `order`
one-pole recursions ``w_j[n] = w_{j-1}[n] + ãw_j[n-1]``, `order` complex multiply–adds
per sample.

The computation runs in the signal's precision: the output element type is `Complex` of
the element type of `x`, whatever precision the filter is stored at.
"""
function gammatone_filt(gf::gammatone_filter, x::AbstractVector{T}) where {T<:AbstractFloat}
    y = Vector{Complex{T}}(undef,length(x))
    return gammatone_filt!(y,gf,x)
end

"""
    gammatone_filt!(y, gf, x)

In-place version of [`gammatone_filt`](@ref): fill the complex vector `y` (same length as
`x`) with the subband signal and return it. Allocation free.
"""
function gammatone_filt!(y::AbstractVector{Complex{T}}, gf::gammatone_filter, x::AbstractVector{<:Real}) where {T<:AbstractFloat}
    length(y) == length(x) || error("Output and input must have the same length")
    ã = convert(Complex{T},gf.coef)
    k = convert(T,gf.norm)
    if(gf.order == 4)
        #Single pass with the four states in registers: the four recursions pipeline
        #across samples, which a stage-by-stage sweep cannot do
        w1 = w2 = w3 = w4 = zero(Complex{T})
        @inbounds for n ∈ eachindex(x,y)
            w1 = muladd(ã,w1,x[n])
            w2 = muladd(ã,w2,w1)
            w3 = muladd(ã,w3,w2)
            w4 = muladd(ã,w4,w3)
            y[n] = k*w4
        end
    else
        #General order: apply the one-pole recursion order times over the signal in place
        @inbounds for n ∈ eachindex(x,y)
            y[n] = x[n]
        end
        for j ∈ 1:gf.order
            w = zero(Complex{T})
            @inbounds for n ∈ eachindex(y)
                w = muladd(ã,w,y[n])
                y[n] = w
            end
        end
        @inbounds for n ∈ eachindex(y)
            y[n] *= k
        end
    end
    return y
end
