#Helper functions: small, general things that are useful across the package and on their own, but
#are not analyses in their own right.
#
#The distinction from spectral_analyses.jl is that nothing here computes a representation - these
#read a number off one, or convert between the ways of describing one. A helper that grows into an
#analysis (returning a spectrum, timefreq or modfreq container, with a comp() form) belongs in the
#file for that analysis instead.

"""
    peak_frequency(x, fs [; <keyword arguments>])
    peak_frequency(s::signal [; <keyword arguments>])
    peak_frequency(sp::spectrum)

Frequency in Hz of the largest component of the magnitude spectrum: the dominant frequency of a
signal. The array and [`signal`](@ref) forms compute a [`magspec`](@ref) first; the
[`spectrum`](@ref) form reads the peak off one already computed, which is the form to use when
the spectrum is wanted for anything else as well.

This is a *convenience over an existing spectrum*, not an estimator. The answer is quantised to
the DFT grid, so it is only as precise as `fs/ndft` — pad with `ndft` for a finer grid, and do not
expect it to resolve two close components or to interpolate between bins. Where a signal's
strongest component is not the one of interest — a tone buried under low-frequency noise, say —
filter first or take the peak over a slice of the spectrum yourself.

### Keyword Arguments
- `ndft` : Number of DFT points [Default is the length of the signal array]
- `wtype` : Window type [Default = "hanning"], see [`window`](@ref)
"""
function peak_frequency(x::AbstractVector{<:AbstractFloat}, fs::Real; ndft::Int = length(x), wtype::String = "hanning")
    mspec,frqs = magspec(comp(), collect(x), fs; ndft, wtype)
    return frqs[argmax(mspec)]
end

peak_frequency(s::signal; ndft::Int = length(s.x), wtype::String = "hanning") =
    peak_frequency(s.x, s.fs; ndft, wtype)

peak_frequency(sp::spectrum) = sp.frqs[argmax(abs.(sp.components))]
