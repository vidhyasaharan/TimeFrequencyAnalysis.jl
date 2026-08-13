#Frequency grid generation and mapping frequencies to grid indices.
#
#Grids are always generated in Float64; the container constructors (types.jl) convert them
#to the precision of the signal they are attached to.

"""
    logfreq_array(; fmin = 10, fmax = 4000, frq_per_octave = 120)

Generate an array of frequencies in Hz, equally spaced in the log domain between `fmin` and
`fmax`, with a resolution of `frq_per_octave` frequencies per octave.
"""
function logfreq_array(;fmin::Real = 10, fmax::Real = 4000, frq_per_octave::Real = 120)
    fmin = convert(Float64,fmin)::Float64
    fmax = convert(Float64,fmax)::Float64
    frq_per_octave = convert(Float64,frq_per_octave)::Float64
    lfmin = log2(fmin)
    lfmax = log2(fmax)
    lfres = 1/frq_per_octave
    lfrq = lfmin:lfres:lfmax
    return exp2.(lfrq)
end

"""
    linfreq_array(; fmin = 0, fmax = 4000, nfrqs = 80)

Generate an array of `nfrqs` equally spaced frequencies in Hz between `fmin` and `fmax`
(both included).
"""
function linfreq_array(;fmin::Real = 0, fmax::Real = 4000, nfrqs::Int = 80)
    fmin = convert(Float64,fmin)::Float64
    fmax = convert(Float64,fmax)::Float64
    fres = (fmax-fmin)/(nfrqs-1)
    frqs = fmin:fres:fmax
    return collect(frqs)
end


## Auditory (ERB) frequency scale

#Glasberg & Moore (1990) ERB constants in the form used by Hohmann (2002), Eq. 13:
#ERB(fc) = ERB_L + fc/ERB_Q
const ERB_L = 24.7  #Hz
const ERB_Q = 9.265 #dimensionless

"""
    erb(fc)

Equivalent rectangular bandwidth (ERB) in Hz of the human auditory filter centred at `fc`
Hz, after Glasberg & Moore (1990): ``ERB(f_c) = 24.7 + f_c/9.265`` (the form used in
Hohmann 2002, Eq. 13). This is the default bandwidth of a [`gammatone_filter`](@ref).
"""
erb(fc::Real) = ERB_L + fc/ERB_Q

"""
    freq2erb(f)

Convert frequency `f` in Hz to its value on the ERB-rate scale, the number of equivalent
rectangular bandwidths below `f`: ``E(f) = 9.265 \\ln(1 + f/(24.7 × 9.265))`` (Hohmann
2002, Eq. 16). Inverse of [`erb2freq`](@ref).
"""
freq2erb(f::Real) = ERB_Q*log(1 + f/(ERB_L*ERB_Q))

"""
    erb2freq(E)

Convert a value `E` on the ERB-rate scale back to frequency in Hz (Hohmann 2002, Eq. 17).
Inverse of [`freq2erb`](@ref).
"""
erb2freq(E::Real) = (ERB_L*ERB_Q)*(exp(E/ERB_Q) - 1)

"""
    erbfreq_array(; fmin = 70, fmax = 6700, base_frq = 1000, filters_per_erb = 1)

Generate an array of frequencies in Hz spaced `1/filters_per_erb` equivalent rectangular
bandwidths apart on the ERB-rate scale (see [`freq2erb`](@ref)): the centre-frequency grid
of a gammatone filterbank in the Hohmann (2002) framework. The defaults reproduce the
~30-channel grid of the reference implementation's demo.

A uniformly spaced grid on the ERB-rate scale is only fully determined once one grid point
is fixed — sliding all points together by less than one spacing would give another equally
valid grid. `base_frq` is that fixed point: the grid consists of every frequency whose
ERB-rate value is `freq2erb(base_frq) + k/filters_per_erb` for integer `k`, clipped to
`[fmin, fmax]`. Choose it as the frequency that matters most in the analysis (a calibration
tone, a fundamental, the 1 kHz reference) and a filter is guaranteed to be centred exactly
there, instead of that frequency falling between two channels. Anchoring at `base_frq`
also means `fmin` and `fmax` only decide how far the grid extends: unlike
[`logfreq_array`](@ref)/[`linfreq_array`](@ref) they are clipping bounds rather than grid
points, and changing them adds or removes edge channels without shifting the interior of
the grid.
"""
function erbfreq_array(;fmin::Real = 70, fmax::Real = 6700, base_frq::Real = 1000, filters_per_erb::Real = 1)
    fmin = convert(Float64,fmin)::Float64
    fmax = convert(Float64,fmax)::Float64
    base_frq = convert(Float64,base_frq)::Float64
    filters_per_erb = convert(Float64,filters_per_erb)::Float64
    (fmin ≤ base_frq ≤ fmax) || error("base_frq must lie within [fmin, fmax]")
    (filters_per_erb > 0) || error("filters_per_erb must be positive")
    ΔE = 1/filters_per_erb
    E0 = freq2erb(base_frq)
    tol = 1e-8 #tolerance in ERB units, so endpoints that lie on the grid survive rounding
    klo = floor(Int, (E0 - freq2erb(fmin))/ΔE + tol)
    khi = floor(Int, (freq2erb(fmax) - E0)/ΔE + tol)
    frqs = [erb2freq(E0 + k*ΔE) for k ∈ -klo:khi]
    frqs[klo+1] = base_frq #pin the anchor exactly (erb2freq∘freq2erb roundtrips only to floating point accuracy)
    return frqs
end


"""
    findclosest(x, data)

Index of the element of `data` closest to `x` (equivalent to `argmin(abs.(data .- x))` but
without allocating). Internal helper behind [`frqindex`](@ref).
"""
function findclosest(x::Real, data::AbstractVector{<:Real})
    T = float(promote_type(typeof(x), eltype(data)))
    d = zero(T)
    mindx::Int = 1
    mmag::T = typemax(T)
    @inbounds for i ∈ eachindex(data)
        d = abs(data[i] - x)
        if(d<mmag)
            mmag = d
            mindx = i
        end
    end
    return mindx
end


"""
    frqindex(f, frqs)

Index into the frequency grid `frqs` of the frequency closest to `f`. When `f` is a vector,
a vector of indices is returned.
"""
frqindex(f::Real, frqs::AbstractVector{<:Real}) = findclosest(f,frqs)

function frqindex(f::AbstractVector{<:Real}, frqs::AbstractVector{<:Real})
    findx = Vector{Int}(undef,length(f))
    for i in eachindex(f)
        findx[i] = frqindex(f[i],frqs)
    end
    return findx
end
