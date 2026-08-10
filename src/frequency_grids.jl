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
