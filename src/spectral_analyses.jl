#Spectral and spectro-temporal analyses: DFT magnitude spectra, STFT spectrograms and
#periodograms evaluated on arbitrary frequency grids.
#
#Each analysis has two forms selected by dispatch (see the comp marker in types.jl):
#the plain form returns a spectrum/timefreq container; the comp() form returns only the
#component array and is the one to use in inner loops.

"""
    dft(x [, <keyword argument>])

Compute the complex DFT coefficients of the array `x` (positive frequencies only, computed
with a real-input FFT). The signal is windowed before transforming.

### Keyword Arguments
- `ndft` : Number of DFT points [Default = `length(x)`]. If `ndft` is less than the length of `x`, the number of DFT points is raised to `length(x)`
- `wtype` : Window type to use [Default = "hanning"], see [`window`](@ref)
"""
function dft(x::AbstractVector{T}; ndft::Int = length(x), wtype::String="hanning") where {T<:AbstractFloat}
    len = length(x)
    win = window(T,len;wtype)
    buflen = max(ndft,len)
    buf = zeros(T, buflen)
    buf[1:len] = x.*win
    return rfft(buf)
end


"""
    magspec(signal [, <keyword arguments>])
    magspec(x, fs = 2π [, <keyword arguments>])

Compute the DFT magnitude spectrum of the [`waveform`](@ref) `signal` (or of the signal in
array `x` with sampling rate `fs`). Returns a [`spectrum`](@ref) object.

    magspec(comp(), x, fs [, <keyword arguments>])

As above, but return the magnitude components and the frequency grid as a tuple of plain
arrays.

    magspec(comp(), x [, <keyword arguments>])

Return only the magnitude components of the DFT of `x`.

### Keyword Arguments
- `ndft` : Number of DFT points [Default is the length of the signal array]
- `wtype` : Window type to use [Default is "hanning"], see [`window`](@ref)
"""
magspec(signal::waveform;ndft::Int = length(signal.x), wtype::String="hanning") = magspec(signal.x, signal.fs; ndft, wtype)

function magspec(x::AbstractVector{<:AbstractFloat},fs::Real=2π; ndft::Int = length(x), wtype::String="hanning")
    mspec,frqs = magspec(comp(), x, fs; ndft, wtype)
    return spectrum(waveform(x,fs), mspec, frqs, "DFT Magnitude Spectrum")
end

function magspec(::comp, x::AbstractVector{<:AbstractFloat},fs::Real; ndft::Int = length(x), wtype::String="hanning")
    ndft = max(ndft, length(x)) #dft() raises the number of DFT points to the signal length; keep the frequency grid consistent with that
    mspec = magspec(comp(), x; ndft, wtype)
    frqs = rfftfreq(ndft, fs)
    return mspec,frqs
end

magspec(::comp, x::AbstractVector{<:AbstractFloat}; ndft::Int = length(x), wtype::String = "hanning") = abs.(dft(x; ndft, wtype))


"""
    specgram(frames::framed_signal [; <keyword arguments>])
    specgram(comp(), frames::framed_signal [; <keyword arguments>])
    specgram([comp(),] x, fs [; <keyword arguments>])

Compute the short-time DFT magnitude spectrogram of the [`framed_signal`](@ref) `frames`
(or of the signal in array `x` with sampling rate `fs`, framed with the durations given by
the keyword arguments). Only the `num_signal_frames` frames that lie entirely within the
signal are analysed. With the [`comp()`](@ref comp) argument the output is the plain matrix
of magnitude components (one column per frame); otherwise it is a [`timefreq`](@ref) object
whose frequency axis matches the rows of the component matrix.

A floor of `eps(T)` is added to every component so that logarithms of the spectrogram are
always well defined.

### Keyword Arguments
- `wtype` : Window type [Default is "hanning"], see [`window`](@ref)
- `ndft` : Number of DFT points [Default is the optimal FFT length that is at least the frame length]
- `frame_dur` : Frame duration in secs, for the signal-array forms [Default = 0.02]
- `frame_shift_dur` : Interval between consecutive frames in secs, for the signal-array forms [Default = 0.01]
"""
function specgram(::comp, sig_frames::framed_signal{T}; ndft::Int = nextfastfft(sig_frames.frame_length), wtype::String = "hanning") where {T<:AbstractFloat}
    len = sig_frames.frame_length
    nfft = max(len,ndft) #Never use fewer DFT points than the frame length
    nframes = sig_frames.num_signal_frames

    win = window(T,len;wtype=wtype)

    buf = zeros(T,nfft) #Zero padded buffer holding one windowed frame at a time
    rfp = plan_rfft(buf) #Precomputed real-input FFT plan, reused for every frame
    nrfft = length(rfp*buf) #Number of (positive) frequency components
    mspec = Matrix{T}(undef,nrfft,nframes)
    for i ∈ 1:nframes
        frame = view_frame(sig_frames,i)
        buf[1:1:len] = win.*frame #Window, then store in the padded buffer
        mspec[:,i] = abs.(rfp*buf)
    end

    #Floor the spectrogram so a log spectrogram is always well defined
    for i ∈ eachindex(mspec)
        mspec[i] += eps(T)
    end
    return mspec
end

function specgram(sig_frames::framed_signal; ndft::Int = nextfastfft(sig_frames.frame_length), wtype::String="hanning")
    mspec = specgram(comp(), sig_frames; ndft, wtype)
    nfft = max(sig_frames.frame_length, ndft) #Number of DFT points used by the computation above
    frqs = rfftfreq(nfft,sig_frames.signal.fs)
    return timefreq(sig_frames,mspec,frqs)
end

function specgram(::comp, x::AbstractVector{<:AbstractFloat}, fs::Real; frame_dur::Real=0.02, frame_shift_dur::Real=0.01, ndft::Int = nextfastfft(time2nsamples(frame_dur,fs)), wtype::String="hanning")
    sig_frames = framed_signal(x,fs,frame_dur,frame_shift_dur)
    return specgram(comp(), sig_frames;ndft, wtype)
end

function specgram(x::Array{<:AbstractFloat},fs::Number;frame_dur::Real=0.02, frame_shift_dur::Real=0.01, ndft::Int = nextfastfft(time2nsamples(frame_dur,fs)), wtype::String="hanning")
    sig_frames = framed_signal(x,fs,frame_dur,frame_shift_dur)
    return specgram(sig_frames; ndft, wtype)
end


## Complex exponential helpers used by the periodogram

"""
    cexp(T, f, fs, N)
    cexp(f, fs, N)

Unit-norm complex negative exponential of length `N` at frequency `f` Hz for sampling rate
`fs`, with element type `Complex{T}` (default `Complex{Float64}`). The periodogram is the
squared magnitude of inner products with these sequences.
"""
function cexp(::Type{T}, f::Real, fs::Real, N::Int) where {T<:AbstractFloat}
    k = T(2π)*(T(f)/T(fs))
    return (1/sqrt(T(N)))*exp.((-im*k).*(1:N))
end

cexp(f::Real,fs::Real,N::Int) = cexp(Float64,f,fs,N)

"""
    cexp_proj_matrix(T, frqs, fs, N)
    cexp_proj_matrix(frqs, fs, N)

Matrix whose row `i` is [`cexp`](@ref)`(T, frqs[i], fs, N)`; multiplying a signal by it
projects the signal onto all the complex exponentials at once.
"""
function cexp_proj_matrix(::Type{T}, frqs::AbstractVector{<:Real},fs::Real,N::Int) where {T<:AbstractFloat}
    nfrqs = length(frqs)
    proj_matrix = zeros(Complex{T},nfrqs,N)
    for i in eachindex(frqs)
        proj_matrix[i,:] = cexp(T,frqs[i],fs,N)
    end
    return proj_matrix
end

cexp_proj_matrix(frqs::AbstractVector{<:Real},fs::Real,N::Int) = cexp_proj_matrix(Float64,frqs,fs,N)


"""
    periodogram([comp(),] x, fs, frqs[; wtype="hanning"])
    periodogram([comp(),] x, fs [; wtype="hanning"[, fmin=10[, fmax=fs/2]]])
    periodogram([comp(),] sig_frames::framed_signal, frqs[; wtype="hanning"])
    periodogram([comp(),] sig_frames::framed_signal[; wtype="hanning"[, fmin=10[, fmax=sig_frames.signal.fs/2]]])

Compute the periodogram of the signal in array `x` with sampling rate `fs` at the
frequencies in `frqs` — or at frequencies equally spaced on the log scale between `fmin` and
`fmax` (see [`logfreq_array`](@ref)) — as the squared magnitude of the inner product between
the windowed signal and a unit-norm complex exponential at each frequency. Unlike the
FFT-based [`magspec`](@ref)/[`specgram`](@ref), the frequency grid is arbitrary.

When the input is a [`framed_signal`](@ref), the periodogram of each frame is computed.
With the [`comp()`](@ref comp) argument the components are returned as a plain
`Vector`/`Matrix`; otherwise a [`spectrum`](@ref)/[`timefreq`](@ref) object is returned.
"""
function periodogram(x::AbstractVector{<:AbstractFloat},fs::Number,frqs::AbstractVector{<:Number};wtype::String="hanning")
    components = periodogram(comp(), x,fs,frqs;wtype)
    signal = waveform(x,fs)
    return spectrum(signal,components,frqs,"Periodogram")
end

function periodogram(x::AbstractVector{<:AbstractFloat},fs::Number;wtype::String="hanning",fmin::Number=10,fmax::Number=fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(x,fs,frqs;wtype=wtype)
end

function periodogram(sig_frames::framed_signal,frqs ;wtype::String="hanning")
    pspec = periodogram(comp(), sig_frames, frqs; wtype)
    return timefreq(sig_frames,pspec,frqs)
end

function periodogram(sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(sig_frames, frqs, wtype = wtype)
end

function periodogram(::comp, x::AbstractVector{T}, fs::Real, frqs::AbstractVector{<:Real}; wtype::String="hanning") where {T<:AbstractFloat}
    len = length(x)
    win = window(T,len;wtype=wtype)
    ip = x.*win
    nfrqs = length(frqs)
    proj = zeros(T,nfrqs)
    for i ∈ eachindex(proj)
        proj[i] = abs2(dot(ip,cexp(T,frqs[i],fs,len)))
    end
    return proj
end

function periodogram(::comp, x::AbstractVector{<:AbstractFloat},fs::Real;wtype::String="hanning",fmin::Real=10,fmax::Real=fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(comp(), x,fs,frqs;wtype=wtype)
end

function periodogram(::comp, sig_frames::framed_signal{T}, frqs ;wtype::String="hanning") where {T<:AbstractFloat}
    len = sig_frames.frame_length
    nframes = sig_frames.num_signal_frames
    fs = sig_frames.signal.fs
    pspec = Matrix{T}(undef,length(frqs),nframes)
    nfrqs = length(frqs)
    win = window(T,len;wtype=wtype)

    #The full projection matrix is precomputed once and reused for every frame; the
    #transpose is collected so that the inner loop walks down contiguous columns
    ce_array = collect(transpose(cexp_proj_matrix(T,frqs,fs,len)))

    for i ∈ 1:nframes
        ip = extract_frame(sig_frames,i)
        ip .*= win
        @views for j in 1:nfrqs
            pspec[j,i] = abs2(dot(ip,ce_array[:,j]))
        end
    end
    return pspec
end

function periodogram(::comp, sig_frames::framed_signal; wtype::String="hanning", fmin=10,fmax=sig_frames.signal.fs/2)
    frqs = logfreq_array(;fmin = fmin,fmax = fmax)
    return periodogram(comp(), sig_frames, frqs, wtype = wtype)
end
