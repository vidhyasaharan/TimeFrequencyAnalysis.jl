#Modulation analysis: the envelope of each filterbank channel, decimated to a rate suited to the
#modulation rates of interest, and the modulation spectrum computed from it.
#
#The chain is signal -> complex subbands -> envelope -> decimation -> modulation transform. Only
#the first two steps know about the filterbank; everything after works on a matrix of envelopes.
#
#The complex subband matrix is never materialised. filt(fb, x) would allocate
#n_channels x n_samples complex numbers - 5.8 GB for an hour of 116-channel sonar - so the
#envelope loop runs one channel at a time through the in-place single-filter kernel, reusing one
#complex buffer, and decimates each row before storing it.

## Envelope rate

#Decimation factor for a requested modulation range. The factor 2.5 rather than 2 keeps
#max_mod_frq clear of the anti-alias filter's transition band: resample places its cutoff at the
#new Nyquist with a transition either side, so asking for fs_env >= 2.5*max_mod_frq leaves the
#requested range inside the flat part rather than on the skirt.
function _decimation_factor(fs::Real, max_mod_frq, fs_env)
    if(fs_env === nothing)
        (max_mod_frq > 0) || error("max_mod_frq must be positive, got $max_mod_frq")
        M = max(1, floor(Int, convert(Float64,fs)/(2.5*convert(Float64,max_mod_frq))))
        return M
    end
    (fs_env > 0) || error("fs_env must be positive, got $fs_env")
    r = convert(Float64,fs)/convert(Float64,fs_env)
    M = round(Int, r)
    if((M < 1) || (abs(r - M) > 1e-9*M))
        lo = fs/max(1,ceil(Int,r))
        hi = fs/max(1,floor(Int,r))
        error("The envelope rate must divide the signal rate by a whole number: fs = $fs and fs_env = $fs_env give a ratio of $r. The nearest usable rates are $lo Hz and $hi Hz")
    end
    return M
end

#Samples at each end of a decimated signal that are still the resampling filter's transient:
#its half-length expressed at the decimated rate, plus one for the rounding.
_resample_margin(::Type{T}, M::Int, rel_bw::Real, attenuation::Real) where {T<:AbstractFloat} =
    ceil(Int, (length(_resample_taps(T, 1//M, rel_bw, attenuation)) - 1)/(2M)) + 1

## Start-up and alignment

"""
    startup_samples(fb::gammatone_filterbank)

Number of samples at the start of a filterbank's output that are still filter start-up rather
than a response to the signal: the slowest channel's [`envelope_delay`](@ref) plus three group
delays, by which point every channel's impulse response has decayed to about 5% of its peak.
[`power_envelope`](@ref) drops this many samples when `trim = true`.
"""
function startup_samples(fb::gammatone_filterbank)
    τ = maximum(group_delay(gf) for gf ∈ fb.filters)
    return maximum(envelope_delay(fb)) + ceil(Int, 3τ)
end

#Turn the align keyword into a gammatone_delay (or nothing). An explicit target shorter than
#default_align leaves the slowest channels clamped at the window edge and silently misaligned,
#which is worth saying out loud - it is the failure mode default_align exists to remove.
_resolve_align(::gammatone_filterbank, ::Nothing) = nothing
_resolve_align(::gammatone_filterbank, d::gammatone_delay) = d

function _resolve_align(fb::gammatone_filterbank, a::Symbol)
    (a === :auto) || error("align must be :auto, a target delay in seconds, a gammatone_delay, or nothing; got :$a")
    return gammatone_delay(fb; delay = :auto)
end

function _resolve_align(fb::gammatone_filterbank, a::Real)
    d = gammatone_delay(fb; delay = a)
    nclamped = count(iszero, d.delays)
    if(nclamped > 0)
        @warn "align = $a s is shorter than default_align(fb) = $(default_align(fb)) s: $nclamped of $(length(d.delays)) channels are too slow to peak inside the search window and stay misaligned"
    end
    return d
end

_envelope_op(envelope::Symbol) = envelope === :power     ? abs2 :
                                 envelope === :amplitude ? abs  :
                                 error("envelope must be :power or :amplitude, got :$envelope")

## Envelopes

"""
    power_envelope([comp(),] s::signal, fb::gammatone_filterbank [; <keyword arguments>])
    power_envelope(comp(), Y::AbstractMatrix{<:Complex}, fs [; <keyword arguments>])

Envelope of every channel of the gammatone filterbank `fb` applied to the signal `s`, decimated
to a sampling rate suited to the modulation rates of interest. Returns a [`timefreq`](@ref)
whose rows are channels and whose columns are time, with `frqs = fb.fcs`; with
[`comp()`](@ref comp) only the component matrix. The second form takes an already-filtered
complex subband matrix (one channel per row) and is `comp()`-only, there being no filterbank or
signal to build a container around.

The envelope is taken from the complex subband signal, so no Hilbert transform or rectify-and-
smooth stage is involved: a gammatone channel is already analytic. Because the filters are
normalised to a peak gain of 2, a tone of amplitude `A` at a channel's centre frequency gives an
amplitude envelope of `A` and a power envelope of `A²`.

### Envelope rate

`max_mod_frq` states the highest modulation rate that must survive, and the decimation factor
follows from it as `M = floor(fs/(2.5·max_mod_frq))`, giving `fs_env = fs/M`. The factor 2.5
rather than 2 keeps `max_mod_frq` clear of the anti-alias filter's transition band. Decimation
is always by a whole number, so the polyphase resampler runs one branch per output sample,
exactly and with linear phase; `fs_env` itself need not be a round number.

Passing `fs_env` overrides `max_mod_frq` and must divide `fs` by a whole number. Anti-aliasing
is done by the resampler (see [`resample`](@ref)), whose filter has unit DC gain, so decimation
leaves the mean of each envelope unchanged.

What decimation discards is modulation above `max_mod_frq`. The opposite limit is set by
the channels themselves: a channel of bandwidth `B` passes both sidebands of a modulation at
`f_m` only while `2f_m ≲ B`, so asking for a modulation rate much beyond half a channel's
bandwidth will find nothing however finely the envelope is sampled.

### Alignment and trimming

Channels differ in delay — at `fc = 77` Hz with a 66 Hz bandwidth the envelope peaks about
7.4 ms late, against well under a millisecond at the top of a bank — so an event appears at a
different time in each row unless that is compensated. `align` is on by default at `:auto`,
which takes its target from [`default_align`](@ref); pass a delay in seconds, a prebuilt
[`gammatone_delay`](@ref), or `nothing` to disable it.

Compensation only shifts and phase-rotates each channel, and the envelope discards phase, so
here it reduces to a per-channel shift. The rows are therefore *sliced* at their own offsets
rather than shifted and zero-filled: no fabricated sample enters the output at all, which is the
same result [`compensate`](@ref)`(...; trim = true)` gives on the complex subbands.

`trim = true` additionally drops [`startup_samples`](@ref)`(fb)` from the head, so no filter
start-up transient survives, and the resampler's own half-length from both ends. The alignment
lead and the start-up allowance are measured from the same origin and do not double-count.

**The time axis keeps the original time origin.** It starts at whatever instant survived
trimming rather than at zero, so the result still refers to the same moments as the input
signal.

It is the time of the *output* sample, the same convention [`gammatone_analysis`](@ref) uses, so
the analysis delay is not removed from it: a modulation event in the signal appears in the
envelope later by the alignment target plus the small difference between the filters' group and
envelope delays (`order/(2πb)` against `(order-1)/(2πb)`, so `1/(2πb)` — about 5 samples at
`fc = 200` Hz with a 93 Hz bandwidth at `fs = 3125`). Without alignment each channel carries its
own delay instead, exactly as it does in the complex subbands. This matters only for reading
absolute event times; `|DFT|²` is shift-invariant, so it does not move anything in a modulation
spectrum.

### Keyword Arguments
- `envelope` : `:power` for `|y|²` or `:amplitude` for `|y|` [Default = `:power`]
- `max_mod_frq` : Highest modulation rate to preserve, in Hz [Default = 125]
- `fs_env` : Envelope sampling rate in Hz, overriding `max_mod_frq` [Default = `nothing`]
- `align` : `:auto`, a target delay in seconds, a [`gammatone_delay`](@ref), or `nothing` [Default = `:auto`]
- `rel_bw` : Relative bandwidth of the anti-alias filter [Default = 1]
- `attenuation` : Stopband attenuation of the anti-alias filter in dB [Default = 60]
- `trim` : Drop the filter start-up and resampler transients [Default = `true`]

# Throws
- `ErrorException` if the filterbank was designed for a different sampling rate, if `envelope`,
  `align`, `max_mod_frq` or `fs_env` is invalid, or if the signal is too short to survive
  trimming.
"""
function power_envelope(::comp, s::signal{T}, fb::gammatone_filterbank;
                        envelope::Symbol = :power, max_mod_frq::Real = 125,
                        fs_env::Union{Nothing,Real} = nothing,
                        align::Union{Nothing,Symbol,Real,gammatone_delay} = :auto,
                        rel_bw::Real = 1, attenuation::Real = 60, trim::Bool = true) where {T<:AbstractFloat}
    check_bank_fs(fb, s)
    op = _envelope_op(envelope)
    M = _decimation_factor(s.fs, max_mod_frq, fs_env)
    d = _resolve_align(fb, align)

    #The alignment lead and the start-up allowance are both measured from the start of the
    #signal, so the head drop is the larger of the two rather than their sum
    lead = (d === nothing) ? 0 : compensation_lead(d)
    head = trim ? max(lead, startup_samples(fb)) : lead

    N = length(s.x)
    nfull = N - head
    (nfull ≥ 1) || error("The signal holds $N samples, of which $head are filter start-up: nothing is left to analyse. Lengthen the signal, or pass trim = false to keep the transient")

    #M = 1 means no decimation at all, so there is no resampling filter and no transient from it
    margin = (trim && M > 1) ? _resample_margin(T, M, rel_bw, attenuation) : 0
    ndec = ceil(Int, nfull/M)
    nout = ndec - 2margin
    (nout ≥ 1) || error("After decimating to $(s.fs/M) Hz only $ndec envelope samples remain, fewer than the $(2margin) the resampling filter's transients occupy. Use a longer signal, a smaller decimation, or trim = false")

    nch = length(fb.filters)
    E = Matrix{T}(undef, nch, nout)
    ybuf = Vector{Complex{T}}(undef, N)
    ebuf = Vector{T}(undef, nfull)
    for k ∈ 1:nch
        filt!(ybuf, fb.filters[k], s.x)
        #Compensation shifts channel k right by Δ and zero-fills the gap; the envelope drops the
        #phase factor, so reading the channel from its own offset gives the same aligned result
        #without inventing a sample
        Δ = (d === nothing) ? 0 : d.delays[k]
        lo = 1 + head - Δ
        @inbounds for i ∈ 1:nfull
            ebuf[i] = op(ybuf[lo + i - 1])
        end
        dec = (M == 1) ? ebuf : resample(ebuf, s.fs, s.fs/M; rel_bw, attenuation)
        @inbounds for j ∈ 1:nout
            E[k,j] = dec[margin + j]
        end
    end
    return E
end

function power_envelope(s::signal{T}, fb::gammatone_filterbank;
                        envelope::Symbol = :power, max_mod_frq::Real = 125,
                        fs_env::Union{Nothing,Real} = nothing,
                        align::Union{Nothing,Symbol,Real,gammatone_delay} = :auto,
                        rel_bw::Real = 1, attenuation::Real = 60, trim::Bool = true) where {T<:AbstractFloat}
    E = power_envelope(comp(), s, fb; envelope, max_mod_frq, fs_env, align, rel_bw, attenuation, trim)
    M = _decimation_factor(s.fs, max_mod_frq, fs_env)
    fsenv = s.fs/M

    #Rebuild the offsets the component form applied, so the time axis reports the instants of the
    #signal these envelope samples came from rather than restarting at zero
    d = _resolve_align(fb, align)
    lead = (d === nothing) ? 0 : compensation_lead(d)
    head = trim ? max(lead, startup_samples(fb)) : lead
    margin = (trim && M > 1) ? _resample_margin(T, M, rel_bw, attenuation) : 0
    t0 = (head + margin*M)/convert(Float64,s.fs)
    time = t0 .+ (0:size(E,2)-1)./convert(Float64,fsenv)

    kind = (envelope === :power) ? "Power" : "Amplitude"
    title = "Gammatone $kind Envelope ($(round(fsenv, digits = 4)) Hz)"
    return timefreq(s, nothing, E, fb.fcs, time, title)
end

function power_envelope(::comp, Y::AbstractMatrix{Complex{T}}, fs::Real;
                        envelope::Symbol = :power, max_mod_frq::Real = 125,
                        fs_env::Union{Nothing,Real} = nothing,
                        align::Union{Nothing,gammatone_delay} = nothing,
                        rel_bw::Real = 1, attenuation::Real = 60, trim::Bool = true) where {T<:AbstractFloat}
    op = _envelope_op(envelope)
    M = _decimation_factor(fs, max_mod_frq, fs_env)

    #No filterbank here, so no start-up allowance can be computed: trim covers the resampler's
    #transients only, and the caller keeps whatever filter start-up the subbands already carry
    lead = (align === nothing) ? 0 : compensation_lead(align)
    (align === nothing) || (size(Y,1) == length(align.delays)) ||
        error("The compensation was built for $(length(align.delays)) channels but the input has $(size(Y,1)) rows")

    N = size(Y,2)
    nfull = N - lead
    (nfull ≥ 1) || error("The subband matrix holds $N samples against an alignment lead of $lead")

    margin = (trim && M > 1) ? _resample_margin(T, M, rel_bw, attenuation) : 0
    ndec = ceil(Int, nfull/M)
    nout = ndec - 2margin
    (nout ≥ 1) || error("After decimating to $(fs/M) Hz only $ndec envelope samples remain, fewer than the $(2margin) the resampling filter's transients occupy")

    nch = size(Y,1)
    E = Matrix{T}(undef, nch, nout)
    ebuf = Vector{T}(undef, nfull)
    for k ∈ 1:nch
        Δ = (align === nothing) ? 0 : align.delays[k]
        lo = 1 + lead - Δ
        @inbounds for i ∈ 1:nfull
            ebuf[i] = op(Y[k, lo + i - 1])
        end
        dec = (M == 1) ? ebuf : resample(ebuf, fs, fs/M; rel_bw, attenuation)
        @inbounds for j ∈ 1:nout
            E[k,j] = dec[margin + j]
        end
    end
    return E
end
