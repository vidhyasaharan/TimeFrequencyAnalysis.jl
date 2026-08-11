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

@testset "waveform recipe" begin
    #A waveform must plot as a single line series of the samples against time in seconds
    #(t = (0:n-1)/fs), with the time axis labelled
    s = waveform(cos.(2pi*440.0 .* (1:800) ./ fs), fs)
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), s)
    @test length(rec) == 1
    t, x = rec[1].args
    @test x == s.x
    @test t == (0:length(s.x)-1) ./ s.fs
    @test rec[1].plotattributes[:seriestype] == :line
    @test rec[1].plotattributes[:xguide] == "Time (secs)"
end

@testset "spectrum recipe" begin
    #A spectrum must plot as a single line series of its components, with x ticks
    #generated from the frequency grid (8 by default, overridable via nticks)
    spec = magspec(signal)
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), spec)
    @test length(rec) == 1
    @test rec[1].args == (spec.components,)
    @test rec[1].plotattributes[:seriestype] == :line
    @test rec[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(spec.frqs, 8)
    #The nticks attribute controls the tick count
    rec4 = RecipesBase.apply_recipe(Dict{Symbol,Any}(:nticks => 4), spec)
    @test rec4[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(spec.frqs, 4)
end

@testset "timefreq recipe" begin
    #A timefreq must plot as a single heatmap series of its component matrix drawn in
    #index space, with tick labels carrying the true time (x) and frequency (y)
    #coordinates
    tf = specgram(framed_signal(signal, 0.02, 0.01))
    rec = RecipesBase.apply_recipe(Dict{Symbol,Any}(), tf)
    @test length(rec) == 1
    @test rec[1].args == (tf.components,)
    @test rec[1].plotattributes[:seriestype] == :heatmap
    @test rec[1].plotattributes[:xticks] == TimeFrequencyAnalysis.generate_ticks(tf.time, 8)
    @test rec[1].plotattributes[:yticks] == TimeFrequencyAnalysis.generate_ticks(tf.frqs, 8)
end
