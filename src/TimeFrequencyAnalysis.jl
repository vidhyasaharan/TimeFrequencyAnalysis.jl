"""
    TimeFrequencyAnalysis

Time-frequency analysis of uniformly sampled signals.

The package provides container types for signals and their spectral representations
([`signal`](@ref), [`framed_signal`](@ref), [`spectrum`](@ref), [`timefreq`](@ref)),
routines to split signals into frames, and spectral / spectro-temporal analyses
([`dft`](@ref), [`magspec`](@ref), [`specgram`](@ref), [`periodogram`](@ref)), together
with supporting utilities (window functions, frequency grids, resampling and synthetic
signal generators). Correlation analyses ([`xcorr`](@ref), [`acorr`](@ref), [`acf`](@ref),
[`nacf`](@ref)) and digital filter responses on arbitrary frequency grids
([`filter_coefs`](@ref), [`filter_resp`](@ref), [`filter_magresp`](@ref)) complete the
classical toolkit. An auditory-model branch adds the complex gammatone filterbank of
Hohmann (2002): ERB-scale frequency grids ([`erbfreq_array`](@ref)), filter and filterbank
design ([`gammatone_filter`](@ref), [`gammatone_filterbank`](@ref)), complex-subband and
cochleagram analyses ([`gammatone_analysis`](@ref), [`gammatone_cochleagram`](@ref)), and
per-channel delay and phase alignment ([`gammatone_delay`](@ref), [`compensate`](@ref)).
Plot recipes for the container types load automatically when the Plots ecosystem
(RecipesBase) is present — a package extension, so none of it is a hard dependency.

All types are parametric in the floating point precision `T<:AbstractFloat` of the signal,
and every routine computes in that precision: a `Float32` signal produces `Float32` results
throughout with no intermediate promotion to `Float64`.

Most analyses come in two flavours, selected by dispatch: called normally they return a
[`spectrum`](@ref)/[`timefreq`](@ref) container carrying the signal, frequency and time
indices; called with [`comp()`](@ref comp) as the first argument they return the plain
component array only.
"""
module TimeFrequencyAnalysis

import DSP                                   #qualified access (DSP.Filters.resample)
import DSP: amp2db, pow2db, filt, filt!      #imported (not just used) so methods for the package's types can be added
using DSP: hamming, hanning, nextfastfft
using FFTW: rfft, plan_rfft, rfftfreq
using LinearAlgebra: dot, mul!
using Random: MersenneTwister

#Signal containers and the component-output marker
export signal, framed_signal, spectrum, timefreq, modfreq, comp

#Framing
export enframe, enframe!, view_frame, extract_frame, number_signal_frames, frame_energy

#Window functions and spectral analyses
export window, dft, magspec, specgram, periodogram, welch_psd

#Element-wise transformations of spectrum/timefreq objects (amp2db and pow2db are DSP.jl
#functions re-exported with added methods; log and log10 methods extend Base)
export amp2db, pow2db

#Correlation sequences and short-time correlation functions
export xcorr, xcorr!, acorr, acorr!, acf, nacf

#Digital filter representation, filtering (filt and filt! are DSP.jl functions re-exported
#with methods for the package's filter types), frequency responses and impulse responses
export filter_coefs, filt, filt!, filter_resp, filter_magresp, filter_impresp

#Frequency grids and frequency indexing (erb and friends are the auditory ERB scale)
export logfreq_array, linfreq_array, frqindex, erb, freq2erb, erb2freq, erbfreq_array

#Gammatone filterbank analysis (Hohmann 2002)
export gammatone_filter, gammatone_impulse_response
export gammatone_filterbank, gammatone_analysis, gammatone_cochleagram
export group_delay, envelope_delay, summed_resp
export gammatone_delay, default_align, compensate, compensate!, compensation_lead

#Sampling helpers
export time2nsamples, resample

#Synthetic signal generators (test and demonstration signals)
export white_noise, ar_process, impulse_train

include("types.jl")             #signal, framed_signal, spectrum, timefreq, comp
include("sampling.jl")          #time2nsamples, resample
include("framing.jl")           #enframe et al., frame_energy, padding helpers
include("windows.jl")           #window
include("frequency_grids.jl")   #logfreq_array, linfreq_array, frqindex, findclosest
include("signal_generators.jl") #white_noise, ar_process, impulse_train
include("spectral_analyses.jl") #dft, magspec, specgram, periodogram, cexp
include("element_ops.jl")       #log/log10/amp2db/pow2db methods for spectrum and timefreq
include("correlations.jl")      #xcorr(!), acorr(!), acf, nacf
include("filters.jl")           #filter_coefs, filter_resp, filter_magresp and helpers
include("gammatone.jl")         #gammatone_filter(bank), gammatone_filt(!), analyses, alignment and helpers
include("plots_support.jl")     #generate_ticks (the recipes live in the RecipesBase extension)

end # module
