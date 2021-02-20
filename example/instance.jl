using ArgParse
using LightGraphs
using MetaGraphs

include("../src/StochasticKep.jl")
using .StochasticKep

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "instance"
        help = "the instance file path"
        arg_type = String
        required = true
    "sample_size"
        help = "SAA sample size"
        arg_type = Int
        required = true
end

parsed_args = parse_args(ARGS, argparse_settings)
instance = parsed_args["instance"]
sample_size = parsed_args["sample_size"]

# ------------------------------------------------- #

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")
wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))

# Original graph
pool = kep_pool(wmd_file, dat_file)

# Matching graph
matching_pool = extract_matching_pool(pool)
to_remove_vertices = [v for v in 1:nv(matching_pool) if length(neighbors(matching_pool, v)) == 0]
for v in reverse(to_remove_vertices)
    rem_vertex!(matching_pool, v)
end

println("*************************************************")
println("Instance : $instance")
println("Relaxed-arc model")
println("-----------------")
println("  Vertices : $(nv(pool))")
println("  Arcs : $(ne(pool))")
println("  First-stage variables : $(ne(pool))")
println("  Second-stage variables : $(sample_size * ne(pool))")
println("Matching model")
println("--------------")
println("  Vertices : $(nv(matching_pool))")
println("  Edges : $(ne(matching_pool))")
println("  First-stage variables : $(ne(matching_pool))")
println("  Second-stage variables : $(sample_size * ne(matching_pool))")
println("*************************************************")


