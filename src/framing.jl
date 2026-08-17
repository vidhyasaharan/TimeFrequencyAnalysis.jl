#Materialising frames from a signal, frame-wise energy, and padding helpers.
#
#A framed_signal (see types.jl) stores only the framing geometry; the functions here produce
#the actual frames, either one at a time (view_frame / extract_frame) or all at once as the
#columns of a matrix (enframe).

"""
    enframe!(A, x)

Fill each column of matrix `A` with a frame of array `x`. The frame size is given by the
number of rows of `A` and the frame shift is chosen so that the final frame ends as close as
possible to the end of `x`.
"""
function enframe!(y::AbstractMatrix{T}, x::AbstractVector{T}) where {T<:AbstractFloat}
    frame_size = size(y,1)
    num_frames = size(y,2)
    len = length(x)
    frame_shift = Int(floor((len - frame_size)/(num_frames-1)))
    for i ∈ 1:num_frames
        sindx = (i-1)*frame_shift + 1
        eindx = sindx + frame_size - 1
        y[:,i] = x[sindx:eindx]
    end
    return y
end


"""
    enframe(f::framed_signal)

Generate a matrix with each column storing one frame of the [`framed_signal`](@ref) `f`,
i.e. materialise every frame that [`view_frame`](@ref) can produce. The matrix holds all
`f.num_frames` frames, so the trailing ones that extend past the end of the signal are zero
padded exactly as `view_frame` returns them.
"""
function enframe(f::framed_signal{T}) where {T<:AbstractFloat}
    y = Matrix{T}(undef,f.frame_length,f.num_frames)
    for i ∈ 1:f.num_frames
        y[:,i] = view_frame(f, i)
    end
    return y
end

"""
    enframe(x, frame_size, frame_shift)

Generate a matrix with each column storing one frame of length `frame_size` samples from the
signal in array `x`, with consecutive frames separated by exactly `frame_shift` samples.
Frames extending past the end of `x` are zero padded, as in [`view_frame`](@ref).
"""
function enframe(x::AbstractVector{T}, frame_size::Int, frame_shift::Int) where {T<:AbstractFloat}
    #Framing a signal read at 1 Hz makes a duration in seconds and a count in samples the
    #same number, so the geometry here is the one framed_signal derives, not a second copy
    return enframe(framed_signal(signal(x,one(T)), frame_size, frame_shift))
end

"""
    enframe(x, fs, frame_dur, frame_shift_dur)

Generate a matrix with each column storing one frame of duration `frame_dur` seconds from
the signal in array `x` with sampling rate `fs`, with consecutive frames separated by
`frame_shift_dur` seconds. Frames extending past the end of `x` are zero padded, as in
[`view_frame`](@ref).
"""
enframe(x::AbstractVector{<:AbstractFloat}, fs::Real, frame_dur::Real, frame_shift_dur::Real) = enframe(framed_signal(signal(x,fs), frame_dur, frame_shift_dur))

"""
    enframe(s::signal, frame_dur, frame_shift_dur)

Generate a matrix with each column storing one frame of duration `frame_dur` seconds from
the [`signal`](@ref) `s`, with consecutive frames separated by `frame_shift_dur`
seconds. Frames extending past the end of the signal are zero padded, as in
[`view_frame`](@ref).
"""
enframe(s::signal, frame_dur::Real, frame_shift_dur::Real) = enframe(framed_signal(s, frame_dur, frame_shift_dur))


"""
    view_frame(x::framed_signal, i)

Create a [view](https://docs.julialang.org/en/v1/manual/performance-tips/#man-performance-views)
of frame number `i` of the [`framed_signal`](@ref) `x`. Frames that extend past the end of
the signal (indices between `num_signal_frames` and `num_frames`) are returned as freshly
allocated zero padded vectors instead of views.
"""
function view_frame(x::framed_signal{T},i::Int) where {T<:AbstractFloat}
    frame_length = x.frame_length
    frame_shift = x.frame_shift
    sindx = (i-1)*frame_shift + 1 #Start index of the desired frame
    eindx = sindx + frame_length - 1 #End index of the desired frame
    if(i<=x.num_signal_frames) #Frame lies entirely within the signal
        frame = @view x.signal.x[sindx:eindx]
    else #Signal ends inside this frame: copy what remains and zero pad to a full frame
        frame = zeros(T,frame_length)
        lindx = length(x.signal.x)
        sig_len = lindx - sindx + 1
        frame[1:sig_len] = @view x.signal.x[sindx:lindx]
    end
    return frame
end


"""
    extract_frame(x::framed_signal, i)

Extract frame number `i` of the [`framed_signal`](@ref) `x` as a new vector (a copying
version of [`view_frame`](@ref)).
"""
extract_frame(x::framed_signal,i::Int) = collect(view_frame(x, i))


"""
    number_signal_frames(x, frame_size, frame_shift)
    number_signal_frames(x, fs, frame_dur, frame_shift_dur)
    number_signal_frames(s::signal, frame_size, frame_shift)
    number_signal_frames(s::signal, frame_dur, frame_shift_dur)

Number of frames that lie entirely within the signal (i.e. without zero padding or
extension) for the array `x` or [`signal`](@ref) `s`, given `frame_size` and `frame_shift`
in samples, or `frame_dur` and `frame_shift_dur` in seconds.
"""
number_signal_frames(x::AbstractVector{<:AbstractFloat}, frame_size::Int, frame_shift::Int) = 1+ Int(floor((length(x)-frame_size)/frame_shift))

number_signal_frames(s::signal, frame_size::Int, frame_shift::Int) = number_signal_frames(s.x, frame_size, frame_shift)

function number_signal_frames(x::AbstractVector{<:AbstractFloat}, fs::Real, frame_dur::Real, frame_shift_dur::Real)
    frame_size = time2nsamples(frame_dur, fs)
    frame_shift = time2nsamples(frame_shift_dur, fs)
    return number_signal_frames(x, frame_size, frame_shift)
end

function number_signal_frames(s::signal, frame_dur::Real, frame_shift_dur::Real)
    frame_size = time2nsamples(frame_dur,s.fs)
    frame_shift = time2nsamples(frame_shift_dur,s.fs)
    return number_signal_frames(s, frame_size, frame_shift)
end


"""
    frame_energy(sig_frames::framed_signal[; normalised = false])

Estimate the energy of each frame of `sig_frames` as the sum of squares of its samples.
With `normalised = true` the energies are scaled so that the maximum is 1. The result
covers all `num_frames` frames (the trailing zero padded ones included).
"""
function frame_energy(sig_frames::framed_signal{T}; normalised = false) where {T<:AbstractFloat}
    numframes = sig_frames.num_frames
    energy = zeros(T,numframes)
    for i=1:numframes
        frame = extract_frame(sig_frames,i)
        energy[i] = sum(abs2,frame)
    end
    if normalised
        return energy/maximum(energy)
    else
        return energy
    end
end


## Padding helpers

"""
    zero_pad(x, pad_len)

Return a copy of vector `x` with `pad_len` zeros prepended and appended.
"""
zero_pad(x::Vector, pad_len::Int) = [zeros(eltype(x),pad_len); x; zeros(eltype(x),pad_len)]

"""
    symmetric_pad(x, pad_len)

Return a copy of vector `x` extended by `pad_len` samples at each end by mirroring the
signal about its end points (the end points themselves are not repeated).
"""
symmetric_pad(x::Vector, pad_len::Int) = [x[pad_len+1:-1:2]; x; x[end-1:-1:end-(pad_len)]]

"""
    unpad_vector(x, pad_len)

Remove `pad_len` samples from each end of vector `x`, undoing [`zero_pad`](@ref) or
[`symmetric_pad`](@ref).
"""
unpad_vector(x::Vector, pad_len::Int) = x[pad_len+1:end-pad_len]
