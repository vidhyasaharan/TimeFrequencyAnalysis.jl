#Window function generation. The underlying window shapes come from DSP.jl; this wrapper
#selects them by name and produces the requested floating point precision.

"""
    window(len [; wtype = "hanning"])
    window(T, len [; wtype = "hanning"])

Generate a window vector of length `len` and type `wtype`, with element type
`T<:AbstractFloat` (default `Float64`). An unrecognised `wtype` falls back to the Hann
window with a warning.

### Implemented window types (`wtype`)
- `"hanning"` : Hann window [default]
- `"hamming"` : Hamming window
- `"rect"` : rectangular window
"""
function window(::Type{T}, flen::Int; wtype::String="hanning") where {T<:AbstractFloat}
    if(wtype=="rect")
        win = ones(T,flen)
    elseif(wtype=="hamming")
        win = convert(Vector{T},hamming(flen))
    elseif(wtype=="hanning")
        win = convert(Vector{T},hanning(flen))
    else
        @warn "Window type \"$wtype\" not recognised - using Hann window"
        win = convert(Vector{T},hanning(flen))
    end
    return win
end

window(flen::Int;wtype::String="hanning") = window(Float64, flen; wtype)
