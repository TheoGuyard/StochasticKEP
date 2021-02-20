using GraphPlot
using LightGraphs
using MetaGraphs

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

instance = "preflib-md/MD-00001-00000020"
wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))

pool = kep_pool(wmd_file, dat_file)
matching_pool = extract_matching_pool(pool)

nodelabel = 1:nv(pool)
gplot(pool, nodelabel=nodelabel, edgelinewidth=0.5)
gplot(matching_pool, nodelabel=nodelabel, edgelinewidth=0.5)