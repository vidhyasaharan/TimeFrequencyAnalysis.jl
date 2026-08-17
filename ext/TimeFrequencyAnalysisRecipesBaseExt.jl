#Plots-ecosystem recipes for the container types, loaded automatically when RecipesBase is
#in the session (Plots loads it). Heatmaps are drawn in index space with generate_ticks
#labels so that nonuniform (e.g. logarithmic) frequency grids render uniformly; no dB or
#other scaling is applied implicitly - pass amp2db(...) etc. explicitly.
module TimeFrequencyAnalysisRecipesBaseExt

using TimeFrequencyAnalysis
using TimeFrequencyAnalysis: generate_ticks
using RecipesBase

#Plot recipe for plotting a spectrogram or other time-frequency representation.
#Everything except the series type is a DEFAULT (-->) rather than forced (:=), so a caller can
#retitle a plot or relabel an axis in the usual way. Forcing the title in particular makes a
#panel of several timefreqs impossible to label, since every panel would carry the same text.
@recipe function f(msp::timefreq; nxticks = 8, nyticks = 8)
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)

    @series begin
        seriestype := :heatmap
        xticks --> xtks
        yticks --> ytks
        xguide --> "Time (sec)"
        yguide --> "Frequency (Hz)"
        if(~isnothing(msp.title))
            title --> msp.title
        end
        msp.components
    end
end


#Plot recipe for plotting a modulation spectrum. Same index-space heatmap as timefreq, since the
#carrier axis is usually ERB or log spaced and would otherwise crush its low channels into the
#bottom of the image; rows are carriers and columns modulation rates, matching the field order.
#Everything except the series type is a DEFAULT (-->) rather than forced (:=), so a caller can
#retitle a plot, relabel an axis or choose another colormap in the usual way. Forcing the title in
#particular makes a panel of several modfreqs impossible to label, since every panel would carry
#the same text.
#
#The default colormap is :ice. A modulation spectrum is a magnitude, so the map has to be
#sequential - one hue, dark to light, with perceived lightness rising monotonically so that equal
#steps in the data look like equal steps on the page. :ice runs near-black to white through blue
#(L* 2 to 98, monotonic), which leaves the noise floor dark and lets only the modulation lines
#carry brightness. A rainbow map would put its brightest colour in the middle of the range, so the
#noise pedestal would compete with the lines for attention - see scripts/modulation-verification.jl
#section 7d/7e, which measures this.
@recipe function f(ms::modfreq; nxticks = 8, nyticks = 8)
    xtks = generate_ticks(ms.mod_frqs,nxticks)
    ytks = generate_ticks(ms.fcs,nyticks)

    @series begin
        seriestype := :heatmap
        seriescolor --> :ice
        xticks --> xtks
        yticks --> ytks
        xguide --> "Modulation Rate (Hz)"
        yguide --> "Carrier Frequency (Hz)"
        if(~isnothing(ms.title))
            title --> ms.title
        end
        ms.components
    end
end


#Plot recipe for plotting spectra. As for timefreq and modfreq, only the series type is forced
#(:=); the title, guide, ticks and the legendless default are defaults (-->) that a caller can
#override - a spectrum drawn over another series may well want its legend back.
@recipe function f(spec::spectrum; nticks = 8)
    frqs = spec.frqs
    xtks = generate_ticks(frqs,nticks)

    legend --> false

    @series begin
        seriestype := :line
        xticks --> xtks
        xguide --> "Frequency(Hz)"
        if(~isnothing(spec.title))
            title --> spec.title
        end
        spec.components
    end
end


#Plot recipe for plotting a time domain signal. Only the series type is forced; the time axis
#label and the legendless default can be overridden by the caller.
@recipe function f(sig::signal)
    x = sig.x
    fs = sig.fs
    t = (0:length(x)-1)/fs

    legend --> false

    @series begin
        seriestype := :line
        xguide --> "Time (secs)"
        t,x
    end
end

end # module
