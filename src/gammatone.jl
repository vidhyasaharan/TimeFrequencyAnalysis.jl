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
    filter_impresp(gf::gammatone_filter, n)
    filter_impresp(fb::gammatone_filterbank, n)

Measured impulse response of the gammatone filter — a unit impulse run through the
filtering kernel for `n` samples — or of a whole filterbank, one channel per row. The
single-filter form matches the closed form of [`gammatone_impulse_response`](@ref) to
machine precision; the bank form is what the alignment measurement
([`gammatone_delay`](@ref)) runs on.
"""
function filter_impresp(gf::gammatone_filter, n::Int)
    (n ≥ 1) || error("The impulse response needs at least one sample")
    e = zeros(Float64,n)
    e[1] = 1.0
    return gammatone_filt(gf,e)
end


"""
    gammatone_impulse_response(gf::gammatone_filter; dur = nothing)

Impulse response of the gammatone filter as a tuple `(t, h)`: sample times `t` in seconds
and the complex response `h`, computed analytically from the exact closed form
``h[n] = k\\binom{n+order-1}{order-1}ã^{\\,n}`` (no filtering involved, so the impulse
response is available directly from the design; [`filter_impresp`](@ref) is the measured
counterpart, which this closed form matches to machine precision). `abs.(h)` is the
response envelope and `real.(h)` the real gammatone impulse response. The duration
defaults to the time the envelope takes to decay to −80 dB relative to its peak; pass
`dur` in seconds to override.
"""
function gammatone_impulse_response(gf::gammatone_filter{T}; dur::Union{Nothing,Real} = nothing) where {T<:AbstractFloat}
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
relative to its peak: the default length used by [`gammatone_impulse_response`](@ref). Found by
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
number of channels); `bw_scale` multiplies whichever is in effect. By default the
bandwidth is realised as the filter's equivalent rectangular bandwidth (the
[`gammatone_filter`](@ref) keyword form); passing `attenuation_db` instead places every
channel's band edges `fc ± B/2` that many dB below its peak (the explicit-attenuation
form, as in pyfilterbank's banks). Note that channel density (`filters_per_erb`) and
bandwidth are independent: one packs more filters, the other widens each of them.

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
        bw::Union{Nothing,Real,AbstractVector{<:Real}} = nothing, bw_scale::Real = 1,
        attenuation_db::Union{Nothing,Real} = nothing)
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
        if(attenuation_db === nothing)
            filters[i] = gammatone_filter(fsT,fcsT[i]; bw = bwi, bw_scale, order)
        else
            #Route B per channel: bandwidth default and scaling as above, edges at
            #fc ± B/2 sitting attenuation_db below the peak
            B = convert(Float64,bw_scale)*((bwi === nothing) ? erb(convert(Float64,fcsT[i])) : convert(Float64,bwi))
            filters[i] = gammatone_filter(fsT,fcsT[i],B,attenuation_db; order)
        end
    end
    return gammatone_filterbank{T}(fsT,fcsT,filters)
end

function gammatone_filterbank(fs::Real; fmin::Real = 70, fmax::Real = min(6700,0.45*fs),
        base_frq::Real = 1000, filters_per_erb::Real = 1, order::Int = 4,
        bw::Union{Nothing,Real,AbstractVector{<:Real}} = nothing, bw_scale::Real = 1,
        attenuation_db::Union{Nothing,Real} = nothing)
    (fmin > 0) || error("fmin must be positive")
    (fmax < fs/2) || error("fmax must lie below fs/2 (got fmax = $fmax with fs/2 = $(fs/2))")
    #The grid is generated in Float64 (like all grids); the bank's precision follows fs
    T = float(typeof(fs))
    frqs = convert(Vector{T},erbfreq_array(;fmin,fmax,base_frq,filters_per_erb))
    return gammatone_filterbank(convert(T,fs),frqs; order, bw, bw_scale, attenuation_db)
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

#Bank method of filter_impresp (documented with the single-filter method above; defined
#here because it needs the gammatone_filterbank type)
function filter_impresp(fb::gammatone_filterbank, n::Int)
    (n ≥ 1) || error("The impulse response needs at least one sample")
    e = zeros(Float64,n)
    e[1] = 1.0
    return gammatone_filt(fb,e)
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
[`gammatone_filterbank`](@ref)`(fs; ...)`. Passing `align` (a target delay in seconds)
folds delay and phase compensation into the analysis: the subbands are aligned with
[`compensate!`](@ref compensate) using [`gammatone_delay`](@ref)`(fb; delay = align)`.
"""
function gammatone_analysis(::comp, s::signal, fb::gammatone_filterbank; align::Union{Nothing,Real} = nothing)
    check_bank_fs(fb,s)
    Y = gammatone_filt(fb,s.x)
    (align === nothing) || compensate!(Y,gammatone_delay(fb;delay = align))
    return Y
end

function gammatone_analysis(s::signal, fb::gammatone_filterbank; align::Union{Nothing,Real} = nothing)
    Y = gammatone_analysis(comp(),s,fb; align)
    t = (0:length(s.x)-1)./s.fs
    title = "Gammatone Filterbank (Hohmann 2002)"
    (align === nothing) || (title *= " (delay compensated)")
    return timefreq(s,Y,fb.fcs,t,title)
end

gammatone_analysis(::comp, s::signal; align::Union{Nothing,Real} = nothing, kwargs...) = gammatone_analysis(comp(),s,gammatone_filterbank(s.fs;kwargs...); align)
gammatone_analysis(s::signal; align::Union{Nothing,Real} = nothing, kwargs...) = gammatone_analysis(s,gammatone_filterbank(s.fs;kwargs...); align)
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
[`gammatone_filterbank`](@ref)`(fs; ...)`. Passing `align` (a target delay in seconds)
delay- and phase-compensates the subbands before the envelope is taken, so the channel
envelopes are time aligned (see [`gammatone_delay`](@ref)).
"""
function gammatone_cochleagram(::comp, s::signal{T}, fb::gammatone_filterbank; align::Union{Nothing,Real} = nothing) where {T<:AbstractFloat}
    Y = gammatone_analysis(comp(),s,fb; align)
    return abs.(Y) .+ eps(T)
end

function gammatone_cochleagram(s::signal, fb::gammatone_filterbank; align::Union{Nothing,Real} = nothing)
    E = gammatone_cochleagram(comp(),s,fb; align)
    t = (0:length(s.x)-1)./s.fs
    title = "Gammatone Cochleagram"
    (align === nothing) || (title *= " (delay compensated)")
    return timefreq(s,E,fb.fcs,t,title)
end

gammatone_cochleagram(::comp, s::signal; align::Union{Nothing,Real} = nothing, kwargs...) = gammatone_cochleagram(comp(),s,gammatone_filterbank(s.fs;kwargs...); align)
gammatone_cochleagram(s::signal; align::Union{Nothing,Real} = nothing, kwargs...) = gammatone_cochleagram(s,gammatone_filterbank(s.fs;kwargs...); align)
gammatone_cochleagram(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_cochleagram(comp(),signal(x,fs);kwargs...)
gammatone_cochleagram(x::AbstractVector{<:AbstractFloat}, fs::Real; kwargs...) = gammatone_cochleagram(signal(x,fs);kwargs...)


## Delay and phase compensation

"""
    gammatone_delay{T<:AbstractFloat}
    gammatone_delay(fb::gammatone_filterbank; delay = 0.004)

Per-channel delay and phase compensation for a filterbank, computed as in `gfb_delay_new`
(Hohmann 2002, §4): the goal is that after compensation every channel's response to an
impulse peaks — envelope and fine structure together — at one common target `delay`
(in seconds, default 4 ms; `nd = round(delay*fs)` samples).

A unit impulse is sent through the bank and each channel's envelope maximum is located
within the window `[0, nd]`; the channel is then delayed by the remaining
`Δ_k = nd - n_peak(k)` samples, and rotated by the unit phase factor ``c_k = i|s_k|/s_k``,
where ``s_k`` is the complex slope of the impulse response at the envelope maximum. At an
interior maximum this rotation makes the channel real and positive at the target instant
(at the envelope peak the slope is purely the fine structure's, ``i`` times the peak
value's phase); channels too slow to peak within the window — normal at the bottom of the
bank — sit at the window edge, get `Δ_k = 0`, and remain partially misaligned, which the
reference design accepts.

Apply with [`compensate`](@ref)/[`compensate!`](@ref compensate), or pass `align` to the
analyses directly.

# Fields
- `fs::T` : sampling rate the compensation was built for, in Hz
- `nd::Int` : target delay in samples
- `delays::Vector{Int}` : per-channel delay ``Δ_k`` in samples
- `phase_factors::Vector{Complex{T}}` : per-channel unit phasors ``c_k``
"""
struct gammatone_delay{T<:AbstractFloat}
    fs::T
    nd::Int
    delays::Vector{Int}
    phase_factors::Vector{Complex{T}}
end

function gammatone_delay(fb::gammatone_filterbank{T}; delay::Real = 0.004) where {T<:AbstractFloat}
    (delay ≥ 0) || error("The alignment target delay must be non-negative")
    nd = Int(round(convert(Float64,delay)*convert(Float64,fb.fs)))
    #Impulse responses over the search window plus two samples, so the slope at a
    #window-edge maximum stays in bounds (the reference implementations do the same)
    H = filter_impresp(fb,nd+3)
    nch = length(fb.filters)
    delays = Vector{Int}(undef,nch)
    phase_factors = Vector{Complex{T}}(undef,nch)
    for k ∈ 1:nch
        m = argmax(abs.(view(H,k,1:nd+1))) #envelope maximum within [0, nd] (m is 1-based)
        delays[k] = nd + 1 - m
        s = (m > 1) ? H[k,m+1] - H[k,m-1] : H[k,2] - H[k,1] #one-sided guard, reachable only for order 1
        phase_factors[k] = convert(Complex{T},im*abs(s)/s)
    end
    return gammatone_delay{T}(fb.fs,nd,delays,phase_factors)
end

"""
    compensate!(Y, d::gammatone_delay)

In-place form of [`compensate`](@ref): apply the per-channel delay and phase compensation
to the complex subband matrix `Y` (one channel per row) and return it.
"""
function compensate!(Y::AbstractMatrix{Complex{T}}, d::gammatone_delay) where {T<:AbstractFloat}
    size(Y,1) == length(d.delays) || error("The compensation was built for $(length(d.delays)) channels but the input has $(size(Y,1)) rows")
    for k ∈ axes(Y,1)
        Δ = d.delays[k]
        c = convert(Complex{T},d.phase_factors[k])
        row = view(Y,k,:)
        N = length(row)
        if(Δ > 0)
            for n ∈ N:-1:(Δ+1)
                row[n] = c*row[n-Δ]
            end
            for n ∈ 1:min(Δ,N)
                row[n] = 0
            end
        else
            for n ∈ 1:N
                row[n] = c*row[n]
            end
        end
    end
    return Y
end

"""
    compensate(Y::AbstractMatrix, d::gammatone_delay)
    compensate(tf::timefreq, d::gammatone_delay)

Apply the per-channel delay and phase compensation `d` to the complex subband
representation returned by [`gammatone_analysis`](@ref):
``ỹ_k[n] = c_k y_k[n - Δ_k]``, zero-filled at the start with the tail truncated, so the
output has the same size as the input. Returns an aligned copy — for a
[`timefreq`](@ref) input, a new `timefreq` with the same signal, frequency and time axes
and the title marked accordingly (see [`compensate!`](@ref) for the in-place matrix
form). Compensation needs the complex analysis output; the cochleagram has already
discarded the phase.
"""
compensate(Y::AbstractMatrix{<:Complex}, d::gammatone_delay) = compensate!(copy(Y),d)

function compensate(tf::timefreq{T,Complex{T}}, d::gammatone_delay) where {T<:AbstractFloat}
    comps = compensate(tf.components,d)
    title = (tf.title === nothing) ? nothing : tf.title*" (delay compensated)"
    return timefreq(tf.signal,tf.frames,comps,tf.frqs,tf.time,title)
end

compensate(::timefreq, ::gammatone_delay) = error("compensate expects the complex gammatone_analysis output; the cochleagram has already discarded the phase")


## Delay and summed-response inspection

"""
    group_delay(gf::gammatone_filter, f)
    group_delay(gf::gammatone_filter)

Group delay of the gammatone filter in **samples** (divide by `fs` for seconds). The
one-argument form returns the closed form at the centre frequency,
``τ_g(f_c) = order⋅λ/(1-λ)`` samples (≈ the analogue ``order/(2πb)`` seconds). The grid
form computes ``τ(f) = -dφ/dθ`` numerically from the unwrapped phase response at each
frequency of `f` in Hz, by central differences — the grid must be dense enough that the
phase moves by less than π between neighbouring points.

Note that the group delay at `fc` is *not* the delay of the envelope maximum: the
asymmetric gammatone envelope peaks earlier (at ≈ ``(order-1)/(2πb)`` seconds against the
group delay's ``order/(2πb)``) — see [`envelope_delay`](@ref).
"""
function group_delay(gf::gammatone_filter)
    λ = convert(Float64,abs(gf.coef))
    return gf.order*λ/(1-λ)
end

function group_delay(gf::gammatone_filter, f::AbstractVector{<:Real})
    length(f) ≥ 2 || error("group_delay needs at least two grid frequencies to differentiate the phase")
    φ = DSP.unwrap(angle.(filter_resp(gf,f)))
    θ = freq2θ(f,gf.fs)
    n = length(f)
    τ = Vector{Float64}(undef,n)
    τ[1] = -(φ[2]-φ[1])/(θ[2]-θ[1])
    τ[n] = -(φ[n]-φ[n-1])/(θ[n]-θ[n-1])
    for i ∈ 2:n-1
        τ[i] = -(φ[i+1]-φ[i-1])/(θ[i+1]-θ[i-1])
    end
    return τ
end

"""
    envelope_delay(gf::gammatone_filter)

Sample index (0-based, so also the delay in samples) of the maximum of the
impulse-response envelope of `gf`, computed from the closed-form impulse response
([`gammatone_impulse_response`](@ref)). This is
the per-channel delay that the alignment stage of the Hohmann framework compensates. The
measured value equals the exact discrete peak position ``⌈(order⋅λ-1)/(1-λ)⌉``, which
sits a sample or two *before* the analogue prediction ``(order-1)/(2πb)`` seconds — the
discrete binomial envelope peaks slightly earlier than its continuous-time counterpart.
"""
function envelope_delay(gf::gammatone_filter)
    _,h = gammatone_impulse_response(gf)
    return argmax(abs.(h)) - 1
end

"""
    summed_resp(fb::gammatone_filterbank, f; delay = nothing)

Complex sum of the channel frequency responses of the filterbank at each frequency of `f`
in Hz: the transfer function of analysing with the bank and summing the raw channel
outputs. Because the channels carry very different delays and phases, the uncompensated
sum interferes and its magnitude ripples deeply.

With `delay` — a [`gammatone_delay`](@ref), or a target delay in seconds from which one is
built — the sum is taken after per-channel compensation,
``Σ_k c_k e^{-iθΔ_k} H_k(e^{iθ})``: the response of analysing, aligning and summing. This
is far flatter; the residual ripple is what the mixer gains of the synthesis stage (v2)
remove.
"""
function summed_resp(fb::gammatone_filterbank, f::AbstractVector{<:Real}; delay::Union{Nothing,Real,gammatone_delay} = nothing)
    if(delay === nothing)
        resp = filter_resp(fb.filters[1],f)
        for i ∈ 2:length(fb.filters)
            resp .+= filter_resp(fb.filters[i],f)
        end
        return resp
    end
    d = (delay isa gammatone_delay) ? delay : gammatone_delay(fb;delay)
    length(d.delays) == length(fb.filters) || error("The delay compensation does not match the number of channels")
    θ = freq2θ(f,fb.fs)
    resp = zeros(ComplexF64,length(f))
    for k ∈ eachindex(fb.filters)
        Hk = filter_resp(fb.filters[k],f)
        ck = convert(ComplexF64,d.phase_factors[k])
        Δ = d.delays[k]
        for i ∈ eachindex(resp)
            resp[i] += ck*cis(-θ[i]*Δ)*Hk[i]
        end
    end
    return resp
end
