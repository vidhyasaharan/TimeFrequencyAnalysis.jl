#Digital filter representation (real or complex coefficients) and frequency responses
#evaluated on arbitrary frequency grids (e.g. the magnitude response of an all-pole AR
#model over a chosen set of analysis frequencies). Transfer functions are evaluated in the
#standard z⁻¹ convention.

"""
    freq2θ(f, fs)

Convert analogue frequency `f` in Hz (scalar or vector) to digital angular frequency
``θ = 2πf/f_s`` in radians per sample, given sampling rate `fs`.
"""
freq2θ(f::Real, fs::Real) = 2pi*(f/fs)

function freq2θ(f::AbstractVector{<:Real}, fs::Real)
    θ = Vector{float(promote_type(eltype(f),typeof(fs)))}(undef,length(f))
    freq2θ!(θ,f,fs)
    return θ
end

"""
    freq2θ!(θ, f, fs)

In-place version of [`freq2θ`](@ref): fill `θ` with the digital angular frequencies of
the analogue frequencies `f` in Hz, given sampling rate `fs`.
"""
function freq2θ!(θ::AbstractVector{<:Real}, f::AbstractVector{<:Real}, fs::Real)
    k = 2pi/fs
    for i ∈ eachindex(θ)
        θ[i] = k*f[i]
    end
    return θ
end

"""
    freq2z(f, fs)

Map analogue frequency `f` in Hz (scalar or vector) to the corresponding point
``z = e^{iθ}`` on the unit circle of the z-plane, given sampling rate `fs`.
"""
freq2z(f::Real,fs::Real) = exp(im*freq2θ(f,fs))

function freq2z(f::AbstractVector{<:Real}, fs::Real)
    z = Vector{Complex{float(promote_type(eltype(f),typeof(fs)))}}(undef,length(f))
    freq2z!(z,f,fs)
    return z
end

"""
    freq2z!(z, f, fs)

In-place version of [`freq2z`](@ref): fill `z` with the unit-circle points corresponding
to the analogue frequencies `f` in Hz, given sampling rate `fs`.
"""
function freq2z!(z::AbstractVector{<:Complex}, f::AbstractVector{<:Real}, fs::Real)
    for i ∈ eachindex(z)
        z[i] = exp(im*freq2θ(f[i],fs))
    end
end

"""
    filter_coefs{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    filter_coefs(num, den)

Coefficients of a rational digital filter: numerator `num` and denominator `den`, constant
terms first — `num[k]` and `den[k]` multiply ``z^{-(k-1)}``. Coefficients may be real or
complex (e.g. the one-pole sections of a complex gammatone filter); `T` is the real
floating point precision and `S` the coefficient element type. The convenience constructor
accepts any real or complex vectors and promotes both to a common element type: all-real
input gives `S = T`, any complex input gives `S = Complex{T}`.
"""
struct filter_coefs{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    num::Vector{S}
    den::Vector{S}
end

function filter_coefs(num::AbstractVector{<:Real}, den::AbstractVector{<:Real})
    T = float(promote_type(eltype(num),eltype(den)))
    return filter_coefs{T,T}(convert(Vector{T},num), convert(Vector{T},den))
end

function filter_coefs(num::AbstractVector{<:Union{Real,Complex}}, den::AbstractVector{<:Union{Real,Complex}})
    T = float(real(promote_type(eltype(num),eltype(den))))
    return filter_coefs{T,Complex{T}}(convert(Vector{Complex{T}},num), convert(Vector{Complex{T}},den))
end


"""
    powers(x, N)

Return the vector `[1, x, x², …, xᴺ]` of the powers of `x` up to `N` (length `N+1`).
"""
function powers(x::Number, N::Int)
    px = Vector{eltype(x)}(undef,N+1)
    px[1] = one(eltype(px))
    for i ∈ 1:N
        px[i+1] = x^i
    end
    return px
end


"""
    H(F, z)
    H(F, f, fs)

Evaluate the rational transfer function described by `F::filter_coefs` at the complex
point `z`, or at analogue frequency `f` in Hz given sampling rate `fs` (i.e. at
``z = e^{i2πf/f_s}``). The numerator and denominator polynomials are evaluated in negative
powers of `z` — the standard ``z^{-1}`` convention — so phase responses and group delays
carry their usual signs (a pure delay ``z^{-1}`` has phase ``-θ``).
"""
function H(F::filter_coefs, z::Complex{<:Real})
    zi = inv(z)
    return evalpoly(zi, F.num)/evalpoly(zi, F.den)
end

H(F::filter_coefs, f::Real, fs::Real) = H(F,freq2z(f,fs))


"""
    filter_resp(F, f, fs)

Complex frequency response of the filter described by `F::filter_coefs` at each analogue
frequency of the vector `f` in Hz, given sampling rate `fs`. See [`H`](@ref) for the
evaluation convention; use [`filter_magresp`](@ref) when only magnitudes are needed.
"""
function filter_resp(F::filter_coefs, f::AbstractVector{<:Real}, fs::Real)
    z = freq2z(f,fs)
    Hz = map(x->H(F,x), z)
    return Hz
end


"""
    Hmag(F, θ)
    Hmag(F, f, fs)

Magnitude of the transfer function described by `F::filter_coefs` at digital angular
frequency `θ` in radians per sample, or at analogue frequency `f` in Hz given sampling
rate `fs`. Equivalent to `abs(H(F, exp(im*θ)))`.
"""
function Hmag(F::filter_coefs{T}, θ::Real) where {T<:AbstractFloat}
    niθ = -im*θ #negative sign: the polynomials are evaluated in powers of z⁻¹ = e^{-iθ}
    CT = Complex{float(promote_type(T,typeof(θ)))}
    Nm = convert(CT,F.num[1])
    Dm = convert(CT,F.den[1])
    Nord = length(F.num) - 1
    Dord = length(F.den) - 1
    if(Nord>0)
        for i ∈ 1:Nord
            Nm += F.num[i+1]*exp(i*niθ)
        end
    end
    if(Dord>0)
        for i ∈ 1:Dord
            Dm += F.den[i+1]*exp(i*niθ)
        end
    end
    return abs(Nm)/abs(Dm)
end

Hmag(F::filter_coefs, f::Real, fs::Real) = Hmag(F, freq2θ(f,fs))


"""
    filter_magresp(F, f, fs)

Magnitude response of the filter described by `F::filter_coefs` at each analogue frequency
of the vector `f` in Hz, given sampling rate `fs`.
"""
filter_magresp(F::filter_coefs, f::AbstractVector{<:Real}, fs::Real) = map(x->Hmag(F,x,fs),f)


"""
    filter_impresp(F, n)

Impulse response of the filter described by `F::filter_coefs`, measured by running a unit
impulse through the filter for `n` samples: `h[k]` is the response at sample `k-1`. The
time-domain complement of [`filter_resp`](@ref)/[`filter_magresp`](@ref) and, like them,
computed in `Float64` regardless of the filter's storage precision. Methods for the
gammatone filter and filterbank types measure through their recursive filtering kernel
instead.
"""
function filter_impresp(F::filter_coefs, n::Int)
    (n ≥ 1) || error("The impulse response needs at least one sample")
    e = zeros(Float64,n)
    e[1] = 1.0
    return filt(F.num,F.den,e)
end
