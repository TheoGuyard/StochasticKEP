using ArgParse, CSV, DataFrames, Dates, JLD2, JuMP, LightGraphs, MetaGraphs, MathOptInterface

include("../src/StochasticKep.jl")
using .StochasticKep

include("typedef.jl")

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "--save"
        help = "save summary as csv file"
        action = :store_true
end

parsed_args = parse_args(ARGS, argparse_settings)

table = DataFrame(
    instance = String[], 
    V = Int[],
    Va  = Int[],
    A = Int[],
    formulation = String[],
    generation = String[],
    budget = Real[],
    K = Int[],
    K_eval = Int[],
    z_sp_perc = Float64[],
    z_sp_real = Float64[],
    z_spr_real = Float64[],
    z_eev_real = Float64[],
    ws = Float64[],
    cycle_1 = Float64[],
    cycle_2 = Float64[],
    cycle_3 = Float64[],
    cycle_4_plus = Float64[],
    t_sp = Float64[],
    t_spr = Float64[],
    gap_sp = Float64[],
    nodes_sp = Float64[],
    simplex_sp = Float64[]
    )

nb_files = length(readdir("benchmarks/results/"))

for (i, file) in enumerate(readdir("benchmarks/results/"))
    println("Reading file $i/$nb_files ...")
    if (split(file, '.')[end] == "jld2") & startswith(file, "MD")
        @load joinpath("benchmarks/results", file) benchmark_instance

        push!(table, [
            split(file, "_")[1],
            nv(benchmark_instance.pool),
            get_prop(benchmark_instance.pool, :nb_altruist),
            ne(benchmark_instance.pool),
            benchmark_instance.formulation,
            benchmark_instance.generation,
            benchmark_instance.budget,
            benchmark_instance.K,
            benchmark_instance.K_eval,
            benchmark_instance.results["sp"]["perc_cost"],
            benchmark_instance.results["sp"]["real_cost"],
            benchmark_instance.results["spr"]["real_cost"],
            benchmark_instance.results["eev"]["real_cost"],
            benchmark_instance.results["ws"]["real_cost"],
            benchmark_instance.results["mean_nb_cycles"]["1"],
            benchmark_instance.results["mean_nb_cycles"]["2"],
            benchmark_instance.results["mean_nb_cycles"]["3"],
            benchmark_instance.results["mean_nb_cycles"]["4+"],
            benchmark_instance.results["sp"]["monitoring"]["solve_time"],
            benchmark_instance.results["spr"]["monitoring"]["solve_time"],
            benchmark_instance.results["sp"]["monitoring"]["relative_gap"],
            benchmark_instance.results["sp"]["monitoring"]["node_count"],
            benchmark_instance.results["sp"]["monitoring"]["simplex_iterations"],
        ])

    end
end

sort!(table,:instance)
println(table)

if parsed_args["save"]
    title = join([
            "benchmarks/saves/",
            Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"),
            ".csv"
        ])
    CSV.write(title, table)
end