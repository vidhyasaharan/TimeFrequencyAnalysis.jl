#Plot recipes from the RecipesBase package extension: the extension must activate
#because this test environment loads RecipesBase, and each recipe must emit the expected
#series (data, series type, ticks, labels). Recipes are exercised through
#RecipesBase.apply_recipe, so no plotting backend is needed.
using RecipesBase

#RecipesBase declares is_key_supported but leaves it to the consuming plotting package
#(Plots defines it); provide it here so recipes with keyword attributes can be applied
#without loading Plots
RecipesBase.is_key_supported(::Symbol) = true

@testset "extension loads" begin
    #Loading RecipesBase alongside TimeFrequencyAnalysis must trigger the extension
    @test Base.get_extension(TimeFrequencyAnalysis, :TimeFrequencyAnalysisRecipesBaseExt) !== nothing
end

@testset "generate_ticks" begin
    #21 uniformly spaced values with 5 ticks: positions must be the evenly spaced indices
    #1, 6, 11, 16, 21 and the labels the coordinates at those indices, rounded to one
    #decimal place
    vals = collect(0.0:0.5:10.0) #21 values
    pos, labels = TimeFrequencyAnalysis.generate_ticks(vals, 5)
    @test pos == [1, 6, 11, 16, 21]
    @test labels == ["0.0", "2.5", "5.0", "7.5", "10.0"]
    #A nonuniform (logarithmic) grid still gets evenly spaced index positions, labelled
    #with the true (log-spaced) coordinates - the point of drawing in index space
    lg = exp10.(range(2.0, 3.0, length = 21))
    lpos, llabels = TimeFrequencyAnalysis.generate_ticks(lg, 3)
    @test lpos == [1, 11, 21]
    @test llabels == string.(round.(lg[lpos], digits = 1))
end

@testset "signal recipe" begin
    #A signal must plot as a single line series of the samples against time in seconds
    #(t = (0:n-1)/fs), with the time axis labelled
    s = signal(cos.(2pi*440.0 .* (1:800) ./ fs), fs)
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), s)
    @test length(rec) == 1
    t, x = rec[1].args
    @test x == s.x
    @test t == (0:length(s.x)-1) ./ s.fs
    @test rec[1].plotattributes[:seriestype] == :line
    @test rec[1].plotattributes[:xguide] == "Time (secs)"

    #The axis label and the legendless default are DEFAULTS, not forced: a caller must be able
    #to relabel the axis and to turn the legend back on (e.g. when overlaying two signals)
    over = RecipesBase.apply_recipe(Dict{Symbol,Any}(:xguide => "t (ms)", :legend => true), s)
    @test over[1].plotattributes[:xguide] == "t (ms)"
    @test over[1].plotattributes[:legend] == true
    #but the series type is not negotiable - a signal is a line
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :scatter), s)[1].plotattributes[:seriestype] == :line
end

@testset "spectrum recipe" begin
    #A spectrum must plot as a single line series of its components, with x ticks
    #generated from the frequency grid (8 by default, overridable via nticks)
    spec = magspec(sig)
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
    @test length(rec) == 1
    @test rec[1].args == (spec.components,)
    @test rec[1].plotattributes[:seriestype] == :line
    @test rec[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(spec.frqs, 8)
    #The nticks attribute controls the tick count
    rec4 = RecipesBase.apply_recipe(Dict{Symbol,Any}(:nticks => 4), spec)
    @test rec4[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(spec.frqs, 4)

    #The title, guide, ticks and legendless default are DEFAULTS, not forced: a caller must be
    #able to retitle a plot, relabel an axis and turn the legend back on
    titled = spectrum(spec.signal, spec.components, spec.frqs, "Magnitude Spectrum")
    over = RecipesBase.apply_recipe(Dict{Symbol,Any}(:title => "panel a", :xguide => "f (Hz)",
                                                     :legend => true), titled)
    @test over[1].plotattributes[:title] == "panel a"
    @test over[1].plotattributes[:xguide] == "f (Hz)"
    @test over[1].plotattributes[:legend] == true
    #an explicit tick spec also wins over the index-space default
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:xticks => ([1,2], ["a","b"])), spec)[1].plotattributes[:xticks] == ([1,2], ["a","b"])
    #but the series type is not negotiable - a spectrum is a line
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :heatmap), spec)[1].plotattributes[:seriestype] == :line
end

@testset "timefreq recipe" begin
    #A timefreq must plot as a single heatmap series of its component matrix drawn in
    #index space, with tick labels carrying the true time (x) and frequency (y)
    #coordinates
    tf = specgram(framed_signal(sig, 0.02, 0.01))
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), tf)
    @test length(rec) == 1
    @test rec[1].args == (tf.components,)
    @test rec[1].plotattributes[:seriestype] == :heatmap
    @test rec[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(tf.time, 8)
    @test rec[1].plotattributes[:yticks] == TimeFrequencyAnalysis.generate_ticks(tf.frqs, 8)

    #A spectrogram is a magnitude, so it defaults to the same sequential colormap the modfreq
    #recipe uses rather than to whatever the plotting package happens to pick - and, being a
    #default, a caller can still choose another
    @test rec[1].plotattributes[:seriescolor] == :ice
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriescolor => :magma), tf)[1].plotattributes[:seriescolor] == :magma

    #The object's title and guides are DEFAULTS, not forced: a caller must be able to retitle a
    #plot and relabel an axis. Without this a panel of several spectrograms cannot be labelled at
    #all, since every panel would carry the object's own title.
    titled = timefreq(tf.frames, tf.components, tf.frqs, "Spectrogram")
    over = RecipesBase.apply_recipe(Dict{Symbol,Any}(:title => "panel a", :xguide => "t (ms)"), titled)
    @test over[1].plotattributes[:title] == "panel a"
    @test over[1].plotattributes[:xguide] == "t (ms)"
    @test over[1].plotattributes[:yguide] == "Frequency (Hz)"   #untouched default survives
    #an explicit tick spec also wins over the index-space default
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:yticks => ([1,2], ["a","b"])), tf)[1].plotattributes[:yticks] == ([1,2], ["a","b"])
    #but the series type is not negotiable - a timefreq is a heatmap
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :line), tf)[1].plotattributes[:seriestype] == :heatmap
end


@testset "modfreq recipe" begin
    #A modfreq must plot as a single heatmap series of its component matrix in index space, with
    #modulation rate on x and carrier frequency on y. The orientation is the part worth pinning:
    #a transposed heatmap still looks like a plausible modulation spectrum, so an asymmetric
    #component matrix is used and the axes checked against the right field.
    ncar, nmod = 5, 9
    C = [10.0k + j for k in 1:ncar, j in 1:nmod]
    fcs = [80.0, 120.0, 180.0, 270.0, 405.0]
    mfr = collect(range(0.0, 40.0; length = nmod))
    ms = modfreq(C, fcs, mfr, 2 .*erb.(fcs), 312.5, 10.0, 11, "Modulation Spectrum")

    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), ms)
    @test length(rec) == 1
    @test rec[1].args == (ms.components,)
    @test rec[1].plotattributes[:seriestype] == :heatmap
    @test rec[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(ms.mod_frqs, 8)
    @test rec[1].plotattributes[:yticks] == TimeFrequencyAnalysis.generate_ticks(ms.fcs, 8)
    @test rec[1].plotattributes[:xguide] == "Modulation Rate (Hz)"
    @test rec[1].plotattributes[:yguide] == "Carrier Frequency (Hz)"
    @test rec[1].plotattributes[:title] == "Modulation Spectrum"

    #The two tick counts are independently controllable, as for timefreq
    rec2 = RecipesBase.apply_recipe(Dict{Symbol,Any}(:nxticks => 4, :nyticks => 3), ms)
    @test rec2[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(ms.mod_frqs, 4)
    @test rec2[1].plotattributes[:yticks] == TimeFrequencyAnalysis.generate_ticks(ms.fcs, 3)

    #The axes are NOT interchangeable: the y ticks must come from fcs, never from mod_frqs
    @test rec[1].plotattributes[:yticks] != TimeFrequencyAnalysis.generate_ticks(ms.mod_frqs, 8)

    #A titleless modfreq emits no title attribute at all, as the timefreq recipe does
    plain = modfreq(C, fcs, mfr, 2 .*erb.(fcs), 312.5, 10.0, 11)
    @test !haskey(RecipesBase.apply_recipe(Dict{Symbol,Any}(), plain)[1].plotattributes, :title)

    #A dB view is still a modfreq and still plots as one
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(), pow2db(ms))[1].plotattributes[:seriestype] == :heatmap

    #The object's title and guides are DEFAULTS, not forced: a caller must be able to retitle a
    #plot and relabel an axis. Without this a panel of several modfreqs cannot be labelled at all,
    #since every panel would carry the object's own title.
    over = RecipesBase.apply_recipe(Dict{Symbol,Any}(:title => "viridis", :xguide => "fₘ (Hz)"), ms)
    @test over[1].plotattributes[:title] == "viridis"
    @test over[1].plotattributes[:xguide] == "fₘ (Hz)"
    @test over[1].plotattributes[:yguide] == "Carrier Frequency (Hz)"   #untouched default survives
    #an explicit tick spec also wins over the index-space default
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:xticks => ([1,2], ["a","b"])), ms)[1].plotattributes[:xticks] == ([1,2], ["a","b"])
    #but the series type is not negotiable - a modfreq is a heatmap
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriestype => :line), ms)[1].plotattributes[:seriestype] == :heatmap

    #The colormap defaults to a sequential, monotonically lightening one rather than to whatever
    #the plotting package happens to use, since a modulation spectrum is a magnitude and a rainbow
    #map would render the noise pedestal as brightly as the lines. It is a default, so a caller
    #can still pick another.
    @test rec[1].plotattributes[:seriescolor] == :ice
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:seriescolor => :magma), ms)[1].plotattributes[:seriescolor] == :magma
    @test RecipesBase.apply_recipe(Dict{Symbol,Any}(:c => :greys), ms)[1].plotattributes[:c] == :greys
end
