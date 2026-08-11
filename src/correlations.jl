#Correlation sequences of discrete-time signals, and short-time correlation functions
#evaluated directly on a waveform (the building blocks of correlation-based pitch
#estimators and of autocorrelation-method AR modelling).

"""
    xcorr(x, h[, z=1])

Compute the cross correlation between `x` and `h`, with the optional `z` indicating the
position of the zero index of the array `h`. The output is of the same length as `x` and
the cross correlation is computed with zero padding.
"""
function xcorr(x::AbstractVector{T}, h::AbstractVector{T}, z::Int=1) where {T}
    y = Vector{T}(undef,length(x))
    xcorr!(y, x, h, z)
    return y
end

"""
    xcorr!(y, x, h[, z=1])

Compute the cross correlation between `x` and `h`, with the optional `z` indicating the
position of the zero index of the array `h`, storing the result in `y`. The number of
sample lags at which the cross correlation is computed is equal to the length of `y` and
the cross correlation is computed with zero padding.
"""
function xcorr!(y::AbstractVector{T}, x::AbstractVector{T}, h::AbstractVector{T}, z::Int=1) where {T}
    padded_x = zeros(T,length(x)+length(h)-1)
    @views padded_x[z:z+length(x)-1] = x
    @views @inbounds for i ∈ eachindex(y)
        y[i] = dot(padded_x[i:i+length(h)-1],h)
    end
end

"""
    acorr!(rxx, x)

Compute the autocorrelation sequence `rxx` of sequence `x`, with the number of sample lags
determined by the length of `rxx`. The autocorrelation value at each lag is normalised by
the length of the correlation window (which reduces at the edges).
"""
function acorr!(rxx::AbstractVector{T}, x::AbstractVector{T}) where {T}
    len = length(x)
    @inbounds @views for i ∈ eachindex(rxx)
        δ = i-1
        N = len - δ
        rxx[i] = (1/N)*dot(x[1:end-δ],x[1+δ:end])
    end
end


"""
    acorr(x, p)

Compute the autocorrelation sequence of `x`, with the number of sample lags given by `p`.
The autocorrelation value at each lag is normalised by the length of the correlation
window (which reduces at the edges).
"""
function acorr(x::AbstractVector{T}, p::Int) where {T}
    rxx = Vector{T}(undef,p)
    acorr!(rxx,x)
    return rxx
end



"""
    acf(sig, i, k; win_size)
    acf(sig, t, lags; win_dur)

Compute the autocorrelation function of the input signal `sig`, provided as a
[`waveform`](@ref) object, at sample `i` (or at time `t` in seconds) at one or more lag
index/indices `k` (or at time lags `lags` in seconds) given a window size of `win_size`
samples (or a window duration `win_dur` in seconds).
"""
acf(s::waveform, i::Int, k::Int; win_size::Int) = dot(s.x[i:i+win_size-k-1], s.x[i+k:i+win_size-1])

function acf(s::waveform, t::Real, lag::Real; win_dur::Real)
    i = time2nsamples(t,s.fs)
    k = time2nsamples(lag,s.fs)
    win_size = time2nsamples(win_dur,s.fs)
    return acf(s,i,k;win_size)
end

acf(s::waveform, i::Int, k::AbstractVector{Int}; win_size::Int) = map(x->acf(s,i,x;win_size),k)

acf(s::waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real) = map(x->acf(s,t,x;win_dur),lags)



"""
    nacf(sig, i, k; win_size, nconst = 0)
    nacf(sig, t, lags; win_dur, nconst = 0)
    nacf(sig; min_lag, max_lag, win_size, win_shift, nconst = 0)

Compute the normalised cross correlation function (with a delayed version of itself) of
the input signal `sig`, provided as a [`waveform`](@ref) object, at sample `i` (or at time
`t` in seconds) at one or more lag index/indices `k` (or at time lags `lags` in seconds)
given a window size of `win_size` samples (or a window duration `win_dur` in seconds).
Both windows are mean-subtracted and the product is scaled by their norms, so values lie
in `[-1, 1]`; a positive `nconst` is added inside the normalising square root to damp the
values of low-energy windows.

The keyword-only form evaluates the function at every lag in `min_lag:max_lag` for
successive windows spaced `win_shift` samples apart, returning a matrix with one column
per window position and one row per lag.
"""
function nacf(s::waveform, i::Int, k::Int; win_size::Int, nconst::Real = 0)
    s1 = s.x[i:i+win_size-1]
    s2 = s.x[i+k:i+k+win_size-1]
    mean_s = sum(s1)/length(s1)
    s1 .-= mean_s
    s2 .-= mean_s
    e1 = sum(abs2,s1)
    e2 = sum(abs2,s2)
    ccf = dot(s1, s2)
    return ccf/(sqrt(nconst + (e1*e2)))
end

function nacf(s::waveform, t::Real, lag::Real; win_dur::Real, nconst::Real = 0.0)
    i = time2nsamples(t, s.fs)
    k = time2nsamples(lag, s.fs)
    win_size = time2nsamples(win_dur, s.fs)
    return nacf(s,i,k;win_size,nconst)
end

nacf(s::waveform, i::Int, k::AbstractVector{Int}; win_size::Int, nconst::Real = 0.0) = map(x->nacf(s,i,x;win_size,nconst),k)

nacf(s::waveform, t::Real, lags::AbstractVector{<:Real}; win_dur::Real, nconst::Real = 0.0) = map(x->nacf(s,t,x;win_dur,nconst),lags)


function nacf(s::waveform{T}; min_lag::Int, max_lag::Int, win_size::Int, win_shift::Int, nconst::Real = 0.0) where {T<:AbstractFloat}
    lags = min_lag:max_lag
    nlags = length(lags)
    nframes = number_signal_frames(s, win_size + max_lag, win_shift)
    cf = Matrix{T}(undef,nlags,nframes)
    for j = 1:nframes
        for i ∈ eachindex(lags)
            sindx = (j-1)*win_shift + 1
            cf[i,j] = nacf(s, sindx, lags[i]; win_size, nconst)
        end
    end
    return cf
end
