#Helpers relating durations in seconds to sample counts, and resampling.

"""
    time2nsamples(dur, fs)

Number of samples spanned by a duration of `dur` seconds at sampling rate `fs` Hz, rounded
to the nearest integer.
"""
time2nsamples(dur::Real, fs::Real) = Int(round(dur*fs))


"""
    resample(s::waveform, fs_new)

Resample the signal `s` to the new sampling rate `fs_new` and return the result as a new
[`waveform`](@ref). Wraps the polyphase `resample` from
[`DSP.jl`](https://docs.juliadsp.org/stable/contents/).
"""
function resample(signal::waveform, fs_new::Number)
    rx = DSP.Filters.resample(signal.x, fs_new/signal.fs)
    return waveform(rx,fs_new)
end
