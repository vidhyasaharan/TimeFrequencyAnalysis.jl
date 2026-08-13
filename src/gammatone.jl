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
    (fc64 > 0.45*fs64) && @warn "Centre frequency $(fc64) Hz is close to the Nyquist frequency; the gammatone design assumes fc ≪ fs/2 and the response skirt will wrap"
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


## Filterbank

"""
    gammatone_filterbank{T<:AbstractFloat}
    gammatone_filterbank(fs; fmin = 70, fmax = min(6700, 0.45fs), base_frq = 1000,
                         filters_per_erb = 1, order = 4, bw = nothing, bw_scale = 1)
    gammatone_filterbank(fs, fcs; order = 4, bw = nothing, bw_scale = 1)

A bank of [`gammatone_filter`](@ref)s sharing the sampling rate `fs`. The centre
frequencies come from [`erbfreq_array`](@ref) — uniformly spaced on the ERB-rate scale
between `fmin` and `fmax` with a channel pinned exactly at `base_frq` (see
[`erbfreq_array`](@ref) for the anchoring logic) — or are given explicitly as the vector
`fcs`. The default `fmax` is the reference implementation's 6700 Hz, capped at `0.45fs`
where the design assumptions end.

Bandwidths follow the single-filter rules per channel: each defaults to its auditory
bandwidth [`erb`](@ref)`(fc)`; `bw` overrides that with one explicit bandwidth for every
channel (scalar — a constant-bandwidth bank) or one per channel (vector matching the
number of channels); `bw_scale` multiplies whichever is in effect. Note that channel
density (`filters_per_erb`) and bandwidth are independent: one packs more filters, the
other widens each of them.

# Fields
- `fs::T` : sampling rate in Hz
- `fcs::Vector{T}` : centre frequency of each channel in Hz
- `filters::Vector{gammatone_filter{T}}` : the channel filters
"""
struct gammatone_filterbank{T<:AbstractFloat}
    fs::T
    fcs::Vector{T}
    filters::Vector{gammatone_filter{T}}
end

function gammatone_filterbank(fs::Real, fcs::AbstractVector{<:Real}; order::Int = 4,
        bw::Union{Nothing,Real,AbstractVector{<:Real}} = nothing, bw_scale::Real = 1)
    isempty(fcs) && error("The filterbank needs at least one centre frequency")
    if(bw isa AbstractVector)
        length(bw) == length(fcs) || error("bw must be a scalar or have one bandwidth per centre frequency")
    end
    T = float(promote_type(typeof(fs),eltype(fcs)))
    fsT = convert(T,fs)
    fcsT = convert(Vector{T},fcs)
    filters = Vector{gammatone_filter{T}}(undef,length(fcsT))
    for i ∈ eachindex(fcsT)
        bwi = (bw isa AbstractVector) ? bw[i] : bw
        filters[i] = gammatone_filter(fsT,fcsT[i]; bw = bwi, bw_scale, order)
    end
    return gammatone_filterbank{T}(fsT,fcsT,filters)
end

function gammatone_filterbank(fs::Real; fmin::Real = 70, fmax::Real = min(6700,0.45*fs),
        base_frq::Real = 1000, filters_per_erb::Real = 1, order::Int = 4,
        bw::Union{Nothing,Real,AbstractVector{<:Real}} = nothing, bw_scale::Real = 1)
    (fmin > 0) || error("fmin must be positive")
    (fmax < fs/2) || error("fmax must lie below fs/2 (got fmax = $fmax with fs/2 = $(fs/2))")
    #The grid is generated in Float64 (like all grids); the bank's precision follows fs
    T = float(typeof(fs))
    frqs = convert(Vector{T},erbfreq_array(;fmin,fmax,base_frq,filters_per_erb))
    return gammatone_filterbank(convert(T,fs),frqs; order, bw, bw_scale)
end

"""
    gammatone_filt(fb::gammatone_filterbank, x)
    gammatone_filt!(Y, fb, x)

Filter the real signal array `x` with every channel of the filterbank, returning the
complex subband matrix with one row per channel: `Y[i,:]` is `x` filtered by
`fb.filters[i]`, and `size(Y) == (length(fb.fcs), length(x))`. The in-place form fills the
preallocated matrix `Y` and returns it. As for the single-filter methods, the computation
runs in the signal's precision.
"""
function gammatone_filt(fb::gammatone_filterbank, x::AbstractVector{T}) where {T<:AbstractFloat}
    Y = Matrix{Complex{T}}(undef,length(fb.filters),length(x))
    return gammatone_filt!(Y,fb,x)
end

function gammatone_filt!(Y::AbstractMatrix{Complex{T}}, fb::gammatone_filterbank, x::AbstractVector{<:Real}) where {T<:AbstractFloat}
    size(Y) == (length(fb.filters),length(x)) || error("Output must be (number of channels) × (signal length)")
    for i ∈ eachindex(fb.filters)
        gammatone_filt!(view(Y,i,:),fb.filters[i],x)
    end
    return Y
end


## Analyses

#The filterbank is designed for one sampling rate; refuse signals at any other
function check_bank_fs(fb::gammatone_filterbank, s::signal)
    isapprox(convert(Float64,fb.fs),convert(Float64,s.fs);rtol=1e-6) || error("The filterbank was designed for fs = $(fb.fs) Hz but the signal has fs = $(s.fs) Hz")
    return nothing
end

"""
    gammatone_analysis([comp(),] s::signal, fb::gammatone_filterbank)
    gammatone_analysis([comp(),] s::signal [; <filterbank keyword arguments>])
    gammatone_analysis([comp(),] x, fs [; <filterbank keyword arguments>])

Complex gammatone filterbank analysis of the [`signal`](@ref) `s` (or of the signal in
array `x` with sampling rate `fs`): every channel of `fb` applied to the signal, giving a
[`timefreq`](@ref) with complex components in which row `i` is the subband signal of the
channel centred at `fcs[i]`, at the signal's full sampling rate (the time axis is
`(0:N-1)/fs` — per-sample resolution, no framing). The magnitude of each row is the band
envelope and its angle the fine structure; see [`gammatone_cochleagram`](@ref) for the
envelope form. With the [`comp()`](@ref comp) argument only the complex component matrix
is returned.

When no filterbank is supplied, one is designed from the keyword arguments via
[`gammatone_filterbank`](@ref)`(fs; ...)`.
"""
function gammatone_analysis(::comp, s::signal, fb::gammatone_filterbank)
    check_bank_fs(fb,s)
    return gammatone_filt(fb,s.x)
end

function gammatone_analysis(s::signal, fb::gammatone_filterbank)
    Y = gammatone_analysis(comp(),s,fb)
    t = (0:length(s.x)-1)./s.fs
    return timefreq(s,Y,fb.fcs,t,"Gammatone Filterbank (Hohmann 2002)")
end

gammatone_analysis(::comp, s::signal; kwargs...) = gammatone_analysis(comp(),s,gammatone_filterbank(s.fs;kwargs...))
gammatone_analysis(s::signal; kwargs...) = gammatone_analysis(s,gammatone_filterbank(s.fs;kwargs...))
gammatone_analysis(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_analysis(comp(),signal(x,fs);kwargs...)
gammatone_analysis(x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_analysis(signal(x,fs);kwargs...)

"""
    gammatone_cochleagram([comp(),] s::signal, fb::gammatone_filterbank)
    gammatone_cochleagram([comp(),] s::signal [; <filterbank keyword arguments>])
    gammatone_cochleagram([comp(),] x, fs [; <filterbank keyword arguments>])

Cochleagram of the [`signal`](@ref) `s` (or of the signal in array `x` with sampling rate
`fs`): the envelope (magnitude) of the [`gammatone_analysis`](@ref), as a real
[`timefreq`](@ref) with one row per channel at the signal's full sampling rate. A floor of
`eps(T)` is added to every component, exactly as in [`specgram`](@ref), so logarithms of
the cochleagram are always well defined (`amp2db(gammatone_cochleagram(s))` plots
directly). With the [`comp()`](@ref comp) argument only the envelope matrix is returned.

When no filterbank is supplied, one is designed from the keyword arguments via
[`gammatone_filterbank`](@ref)`(fs; ...)`.
"""
function gammatone_cochleagram(::comp, s::signal{T}, fb::gammatone_filterbank) where {T<:AbstractFloat}
    Y = gammatone_analysis(comp(),s,fb)
    return abs.(Y) .+ eps(T)
end

function gammatone_cochleagram(s::signal, fb::gammatone_filterbank)
    E = gammatone_cochleagram(comp(),s,fb)
    t = (0:length(s.x)-1)./s.fs
    return timefreq(s,E,fb.fcs,t,"Gammatone Cochleagram")
end

gammatone_cochleagram(::comp, s::signal; kwargs...) = gammatone_cochleagram(comp(),s,gammatone_filterbank(s.fs;kwargs...))
gammatone_cochleagram(s::signal; kwargs...) = gammatone_cochleagram(s,gammatone_filterbank(s.fs;kwargs...))
gammatone_cochleagram(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_cochleagram(comp(),signal(x,fs);kwargs...)
gammatone_cochleagram(x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_cochleagram(signal(x,fs);kwargs...)
