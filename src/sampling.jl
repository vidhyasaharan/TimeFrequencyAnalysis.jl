#Helpers relating durations in seconds to sample counts, and resampling.

"""
    time2nsamples(dur, fs)

Number of samples spanned by a duration of `dur` seconds at sampling rate `fs` Hz, rounded
to the nearest integer.
"""
time2nsamples(dur::Real, fs::Real) = Int(round(dur*fs))


"""
    resample(s::signal, fs_new)

Resample the signal `s` to the new sampling rate `fs_new` and return the result as a new
[`signal`](@ref). Wraps the polyphase `resample` from
[`DSP.jl`](https://docs.juliadsp.org/stable/contents/).
"""
function resample(s::signal, fs_new::Number)
    rx = DSP.Filters.resample(s.x, fs_new/s.fs)
    return signal(rx,fs_new)
end
