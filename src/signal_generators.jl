#Synthetic signal generators, mainly for tests, demonstrations and quick experiments.

"""
    white_noise([T,] len)
    white_noise(dur, fs)

Generate `len` samples (or `dur` seconds at sampling rate `fs`) of unit-variance Gaussian
white noise with element type `T<:AbstractFloat` (default `Float64`). Each call draws from a
freshly seeded random number generator.
"""
white_noise(::Type{T}, len::Int) where {T<:AbstractFloat} = randn(MersenneTwister(), T, len)
white_noise(len::Int) = white_noise(Float64, len)
white_noise(dur::Real, fs::Real) = white_noise(time2nsamples(dur,fs))

"""
    ar_process(a, len)
    ar_process(a, dur, fs)

Generate `len` samples (or `dur` seconds at sampling rate `fs`) of an autoregressive process
by filtering white noise with the all-pole filter `1/A(z)` whose denominator coefficients
are given by the vector `a` (ordered as `[1, a₁, a₂, ...]`).
"""
ar_process(a::Vector, len::Int) = filt([1],a, white_noise(len))
ar_process(a::Vector, dur::Real, fs::Real) = ar_process(a,time2nsamples(dur,fs))

"""
    impulse_train([T,] period, len)
    impulse_train(f₀, dur, fs)

Generate `len` samples of a unit impulse train with `period` samples between impulses
(the first impulse is at sample 1), with element type `T<:AbstractFloat` (default
`Float64`). Alternatively specify an impulse rate `f₀` in Hz and a duration `dur` in seconds
at sampling rate `fs`.
"""
function impulse_train(::Type{T}, period::Int, len::Int) where {T<:AbstractFloat}
    x = zeros(T,len)
    x[1:period:end] .= 1
    return x
end

impulse_train(period::Int, len::Int) = impulse_train(Float64, period, len)

impulse_train(f₀::Real, dur::Real, fs::Real) = impulse_train(Int(round(fs/f₀)), Int(round(dur*fs)))
