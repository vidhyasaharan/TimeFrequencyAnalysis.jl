#Element-wise transformations of spectrum and timefreq objects.
#
#Each method applies a scalar function to every stored component and returns a new object of
#the same kind with all other fields (signal, framing, frequency, time and modulation indices,
#title) unchanged. log and log10 extend Base; amp2db and pow2db extend DSP.jl.

"""
    log(sp::spectrum) / log(tf::timefreq) / log(mf::modfreq)

Elementwise natural logarithm of the stored components, returned as a new
[`spectrum`](@ref)/[`timefreq`](@ref)/[`modfreq`](@ref) with all other fields unchanged.
"""
Base.log(sp::spectrum) = spectrum(sp.signal, map(log, sp.components), sp.frqs, sp.title)
Base.log(tf::timefreq) = timefreq(tf.signal, tf.frames, map(log, tf.components), tf.frqs, tf.time, tf.title)
Base.log(mf::modfreq) = modfreq(map(log, mf.components), mf.fcs, mf.mod_frqs, mf.bws, mf.fs_env, mf.frame_dur, mf.nframes, mf.title)

"""
    log10(sp::spectrum) / log10(tf::timefreq) / log10(mf::modfreq)

Elementwise base-10 logarithm of the stored components, returned as a new
[`spectrum`](@ref)/[`timefreq`](@ref)/[`modfreq`](@ref) with all other fields unchanged.
"""
Base.log10(sp::spectrum) = spectrum(sp.signal, map(log10, sp.components), sp.frqs, sp.title)
Base.log10(tf::timefreq) = timefreq(tf.signal, tf.frames, map(log10, tf.components), tf.frqs, tf.time, tf.title)
Base.log10(mf::modfreq) = modfreq(map(log10, mf.components), mf.fcs, mf.mod_frqs, mf.bws, mf.fs_env, mf.frame_dur, mf.nframes, mf.title)

"""
    amp2db(sp::spectrum) / amp2db(tf::timefreq) / amp2db(mf::modfreq)

Convert amplitude components to decibels (`20log10`) elementwise, returned as a new
[`spectrum`](@ref)/[`timefreq`](@ref)/[`modfreq`](@ref) with all other fields unchanged. Extends
`amp2db` from DSP.jl.
"""
amp2db(sp::spectrum) = spectrum(sp.signal, map(amp2db, sp.components), sp.frqs, sp.title)
amp2db(tf::timefreq) = timefreq(tf.signal, tf.frames, map(amp2db, tf.components), tf.frqs, tf.time, tf.title)
amp2db(mf::modfreq) = modfreq(map(amp2db, mf.components), mf.fcs, mf.mod_frqs, mf.bws, mf.fs_env, mf.frame_dur, mf.nframes, mf.title)

"""
    pow2db(sp::spectrum) / pow2db(tf::timefreq) / pow2db(mf::modfreq)

Convert power components to decibels (`10log10`) elementwise, returned as a new
[`spectrum`](@ref)/[`timefreq`](@ref)/[`modfreq`](@ref) with all other fields unchanged. Extends
`pow2db` from DSP.jl.
"""
pow2db(sp::spectrum) = spectrum(sp.signal, map(pow2db, sp.components), sp.frqs, sp.title)
pow2db(tf::timefreq) = timefreq(tf.signal, tf.frames, map(pow2db, tf.components), tf.frqs, tf.time, tf.title)
pow2db(mf::modfreq) = modfreq(map(pow2db, mf.components), mf.fcs, mf.mod_frqs, mf.bws, mf.fs_env, mf.frame_dur, mf.nframes, mf.title)
