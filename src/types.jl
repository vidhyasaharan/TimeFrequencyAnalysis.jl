#Container types shared by all analyses in the package.
#
#Every container is parametric in T<:AbstractFloat, the floating point precision of the
#signal it was built from. Constructors convert their numeric inputs to that precision, so
#analyses never mix precisions internally.

"""
    signal(x, fs)

Store the signal in array `x` with sampling rate `fs` (in Hz) as a `signal` object.

The element type of `x` (any `AbstractFloat`) sets the numeric precision used by all
subsequent processing of the signal. If `x` is a matrix it is assumed to hold one channel
per row or per column (whichever makes the channels longer) and the first channel is stored;
a warning is emitted in that case.

# Fields
- `x::Vector{T}` : signal samples
- `fs::T` : sampling rate in Hz
"""
struct signal{T<:AbstractFloat}
    x::Vector{T}
    fs::T
end

function signal(x::AbstractArray{T}, fs::Number) where {T<:AbstractFloat}
    if(ndims(x)>2)
        error("Signal has more than 2 dimensions, interpretation is not known")
    end
    if(ndims(x)==2)
        if((size(x,1)==1)||(size(x,2)==1))
            samples = x[:]
        else
            @warn "Input is a matrix; assuming one channel per $(size(x,1)>size(x,2) ? "column" : "row") and taking the first channel"
            if(size(x,1)>size(x,2))
                samples = x[:,1]
            else
                samples = x[1,:]
            end
        end
    else
        samples = x
    end
    samples = convert(Vector{T},samples)
    return signal(samples,convert(T,fs))
end


"""
    framed_signal(sig::signal, win_dur=0.02, win_shift=0.01)
    framed_signal(x, fs, win_dur=0.02, win_shift=0.01)

Describe a division of the signal in the [`signal`](@ref) `sig` (or in array `x` with
sampling rate `fs`) into frames of duration `win_dur` seconds whose start points are
`win_shift` seconds apart. No copies of the signal are made: the object stores the signal
plus the framing geometry, and frames are materialised on demand by [`view_frame`](@ref),
[`extract_frame`](@ref) or [`enframe`](@ref).

# Fields
- `signal::signal{T}` : the underlying signal
- `frame_length::Int` : frame length in samples
- `frame_shift::Int` : interval between the start of consecutive frames in samples
- `num_signal_frames::Int` : number of frames that lie entirely within the signal
- `num_frames::Int` : number of frames needed to cover every sample, where the final frames
  may extend past the end of the signal and are zero padded when materialised
"""
struct framed_signal{T<:AbstractFloat}
    signal::signal{T}
    frame_length::Int
    frame_shift::Int
    num_signal_frames::Int
    num_frames::Int
end

function framed_signal(s::signal,win_dur::Number=0.02,win_shift::Number=0.01)
    fs = s.fs
    x = s.x
    frame_length = Int(round(win_dur*fs))
    frame_shift = Int(round(win_shift*fs))
    num_signal_frames = Int(floor(1 + (length(x)-frame_length)/frame_shift))
    num_frames = Int(ceil(length(x)/frame_shift))
    return framed_signal(s,frame_length,frame_shift,num_signal_frames,num_frames)
end

function framed_signal(x::Array{<:AbstractFloat},fs::Number,win_dur::Number=0.02,win_shift::Number=0.01)
    s = signal(x,fs)
    return framed_signal(s,win_dur,win_shift)
end


"""
    spectrum(s, components, frqs, title)
    spectrum(s, components, frqs)

Hold the spectral `components` (real or complex, matching the precision of `s`) of the
[`signal`](@ref) `s`, where component `i` corresponds to the frequency `frqs[i]` in
Hz. `title` is free text describing the analysis that produced the spectrum (defaults to
`""`) and is used, for example, by plot recipes.

# Fields
- `signal::signal{T}` : the analysed signal
- `components::Vector{S}` : spectral components, `S` is `T` or `Complex{T}`
- `frqs::Vector{T}` : frequency of each component in Hz
- `title::AbstractString` : description of the analysis
"""
struct spectrum{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    signal::signal{T}
    components::Vector{S}
    frqs::Vector{T}
    title::AbstractString
end

function spectrum(s::signal{T}, components::AbstractVector{S}, frqs::AbstractVector{<:Real}, title::AbstractString) where {T<:AbstractFloat, S<:Union{T,Complex{T}}}
    return spectrum(s, convert(Vector{S},components), convert(Vector{T},frqs), title)
end

spectrum(s::signal, components::AbstractVector, frqs::AbstractVector{<:Real}) = spectrum(s,components,frqs,"")


"""
    timefreq(s, frames, components, frqs, time, title)
    timefreq(s, components, frqs, time, title)
    timefreq(frames, components, frqs, time, title)
    timefreq(frames, components, frqs, title)
    timefreq(frames, components, frqs)

Hold a time-frequency (spectro-temporal) representation of a signal: a matrix of
`components` in which each column is a spectral snapshot at one time instant and each row
tracks one frequency over time.

The representation is anchored to the [`signal`](@ref) `s` it was computed from and,
when it was obtained frame-by-frame, to the [`framed_signal`](@ref) `frames` describing the
framing (`frames` is `nothing` otherwise). When `time` is omitted, frame centre times are
computed from the framing geometry.

# Fields
- `signal::signal{T}` : the analysed signal
- `frames::Union{framed_signal{T}, Nothing}` : framing used, if any
- `components::Matrix{S}` : time-frequency components, `S` is `T` or `Complex{T}`;
  `size(components) == (length(frqs), number of time instants)`
- `frqs::Vector{T}` : frequency of each row in Hz
- `time::Union{Vector{T}, Vector{Vector{T}}}` : time of each column in seconds, or one time
  vector per row when rows are not sampled at common instants
- `title::Union{AbstractString, Nothing}` : description of the analysis
"""
struct timefreq{T<:AbstractFloat, S<:Union{T,Complex{T}}}
    signal::signal{T}
    frames::Union{framed_signal{T}, Nothing}
    components::Matrix{S}
    frqs::Vector{T}
    time::Union{Vector{T}, Vector{Vector{T}}}
    title::Union{AbstractString, Nothing}
end

#Convert time indices (flat vector, or vector of vectors) to the precision of the signal
_convert_time(::Type{T}, time::AbstractVector{<:Real}) where {T<:AbstractFloat} = convert(Vector{T},time)
_convert_time(::Type{T}, time::AbstractVector{<:AbstractVector{<:Real}}) where {T<:AbstractFloat} = [convert(Vector{T},t) for t in time]

#Central converting constructor: all other constructors funnel into this one, which converts
#every field to the precision of the signal and then calls the default constructor
function timefreq(s::signal{T},
        frames::Union{framed_signal{T}, Nothing},
        components::AbstractMatrix{S},
        frqs::AbstractVector{<:Real},
        time::Union{AbstractVector{<:Real}, AbstractVector{<:AbstractVector{<:Real}}},
        title::Union{AbstractString, Nothing}) where {T<:AbstractFloat, S<:Union{T,Complex{T}}}
    return timefreq(s,frames,convert(Matrix{S},components),convert(Vector{T},frqs),_convert_time(T,time),title)
end

timefreq(frames::framed_signal,
        components::AbstractMatrix,
        frqs::AbstractVector{<:Real},
        time::Union{AbstractVector{<:Real}, AbstractVector{<:AbstractVector{<:Real}}},
        title::Union{AbstractString, Nothing}) = timefreq(frames.signal,frames,components,frqs,time,title)

timefreq(s::signal,
        components::AbstractMatrix,
        frqs::AbstractVector{<:Real},
        time::Union{AbstractVector{<:Real}, AbstractVector{<:AbstractVector{<:Real}}},
        title::Union{AbstractString, Nothing}) = timefreq(s,nothing,components,frqs,time,title)

#When no time indices are given, place each column at the centre time of its frame
function timefreq(frames::framed_signal, components::AbstractMatrix, frqs::AbstractVector{<:Real}, title::Union{AbstractString, Nothing})
    frame_shift = frames.frame_shift
    fs = frames.signal.fs
    frame_length = frames.frame_length
    num_frames = size(components,2)

    time_shift = frame_shift/fs
    start_time = frame_length/(2*fs)
    t = collect(range(start_time, step = time_shift, length = num_frames))
    return timefreq(frames.signal,frames,components,frqs,t,title)
end

timefreq(frames::framed_signal, components::AbstractMatrix, frqs::AbstractVector{<:Real}) = timefreq(frames,components,frqs,nothing)


"""
    comp

Dispatch marker requesting plain component output. Analyses that normally return a
[`spectrum`](@ref) or [`timefreq`](@ref) object return only the underlying component
`Vector`/`Matrix` when `comp()` is passed as their first argument:

```julia
tf = specgram(frames)           #timefreq object
S  = specgram(comp(), frames)   #the same components as a plain Matrix
```
"""
struct comp end
