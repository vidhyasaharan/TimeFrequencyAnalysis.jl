#Plotting support that is ordinary code. The Plots-ecosystem recipes for
#waveform/spectrum/timefreq live in the package extension
#ext/TimeFrequencyAnalysisRecipesBaseExt.jl and activate when RecipesBase is loaded
#(in practice, when Plots is).

"""
    generate_ticks(label_values, nticks)

Build a Plots-style tick tuple `(positions, labels)` with `nticks` ticks for an axis whose
data is indexed `1:length(label_values)` but represents the coordinates stored in
`label_values` (a time or frequency array): positions are evenly spaced indices, and each
label is the corresponding coordinate rounded to one decimal place. This lets a heatmap
drawn in index space — which renders nonuniform (e.g. logarithmic) grids uniformly —
display true coordinates on its axes.

Unexported; used by the RecipesBase package extension and by downstream packages' recipes.
"""
function generate_ticks(label_values::AbstractVector{<:Real}, nticks::Int)
    indx = Int.(round.(range(1, length(label_values), length = nticks)))
    tks = (indx,string.(round.(label_values[indx],digits=1)))
    return tks
end
