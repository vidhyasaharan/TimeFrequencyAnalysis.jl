#Plots-ecosystem recipes for the container types, loaded automatically when RecipesBase is
#in the session (Plots loads it). Heatmaps are drawn in index space with generate_ticks
#labels so that nonuniform (e.g. logarithmic) frequency grids render uniformly; no dB or
#other scaling is applied implicitly - pass amp2db(...) etc. explicitly.
module TimeFrequencyAnalysisRecipesBaseExt

using TimeFrequencyAnalysis
using TimeFrequencyAnalysis: generate_ticks
using RecipesBase

#Plot recipe for plotting a spectrogram or other time-frequency representation
@recipe function f(msp::timefreq; nxticks = 8, nyticks = 8)
    frqs = msp.frqs
    time = msp.time
    xtks = generate_ticks(time,nxticks)
    ytks = generate_ticks(frqs,nyticks)

    @series begin
        seriestype := :heatmap
        xticks := xtks
        yticks := ytks
        xguide := "Time (sec)"
        yguide := "Frequency (Hz)"
        if(~isnothing(msp.title))
            title := msp.title
        end
        msp.components
    end
end


#Plot recipe for plotting spectra
@recipe function f(spec::spectrum; nticks = 8)
    frqs = spec.frqs
    xtks = generate_ticks(frqs,nticks)

    legend := false

    @series begin
        seriestype := :line
        xticks := xtks
        xguide := "Frequency(Hz)"
        if(~isnothing(spec.title))
            title := spec.title
        end
        spec.components
    end
end


#Plot recipe for plotting a time domain signal
@recipe function f(sig::signal)
    x = sig.x
    fs = sig.fs
    t = (0:length(x)-1)/fs

    legend := false

    @series begin
        seriestype := :line
        xguide := "Time (secs)"
        t,x
    end
end

end # module
