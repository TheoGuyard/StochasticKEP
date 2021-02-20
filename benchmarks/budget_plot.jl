using ArgParse, CSV, DataFrames, JLD2, MetaGraphs, Plots

include("../src/StochasticKep.jl")
using .StochasticKep

include("typedef.jl")

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "file"
        help = "the budget benchmark file"
        arg_type = String
        required = true
end

parsed_args = parse_args(ARGS, argparse_settings)
file = parsed_args["file"]

instance = split(file, "_")[2]

println("Reading budget instance $instance ...")

@load joinpath("benchmarks/results", file) budget_instance
bi = budget_instance

p1 = plot(
    bi.budget_factors,
    [bi.z_sp_perc, bi.z_sp_real, bi.eev, bi.ws],
    label = ["perc cost" "real cost" "EEV" "WS"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    xlabel = "Budget factor",
    size = (600, 400)
)

p2 = plot(
    bi.budget_factors,
    bi.t_sp_perc,
    label = "",
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :topright,
    xlabel = "Budget factor",
    ylabel = "Solution time (seconds)",
    size = (600, 400)
)

title = join([instance, bi.generation], " | ")

plot(p1, p2, 
    layout=(1,2), 
    title=title, 
    size = (1200, 400), 
    bottom_margin = 5Plots.mm
    )