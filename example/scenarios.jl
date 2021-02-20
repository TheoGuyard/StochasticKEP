using ArgParse
using DataFrames
using LightGraphs
using Plots
using Statistics

include("../src/StochasticKep.jl")
using .StochasticKep

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "instance"
        help = "the instance file path"
        arg_type = String
        required = true
    "formulation"
        help = "stochastic model formulation"
        arg_type = String
        required = true
    "generation"
        help = "failure generation type"
        arg_type = String
        required = true
    "budget"
        help = "HLA-crossmatch test budget"
        arg_type = Int
        required = true
    "min_sample_size"
        help = "minimum SAA sample size"
        arg_type = Int
        required = true
    "max_sample_size"
        help = "maximum SAA sample size"
        arg_type = Int
        required = true
    "step"
        help = "step in SAA sample size"
        arg_type = Int
        required = true
    "eval_sample_size"
        help = "number of SAA samples in the evaluation of a decision"
        arg_type = Int
        required = true
    "repeats"
        help = "number of repeats per sample size"
        arg_type = Int
        required = false
        default = 1
    "maxtime"
        help = "maximum solving time"
        arg_type = Int
        required = false
        default = 60
end

parsed_args = parse_args(ARGS, argparse_settings)
instance = parsed_args["instance"]
formulation = parsed_args["formulation"]
generation = parsed_args["generation"]
budget = parsed_args["budget"]
min_sample_size = parsed_args["min_sample_size"]
max_sample_size = parsed_args["max_sample_size"]
step = parsed_args["step"]
eval_sample_size = parsed_args["eval_sample_size"]
repeats = parsed_args["repeats"]
maxtime = parsed_args["maxtime"]

budget > 0 || throw(ArgumentError("Excpected budget > 0."))
max_sample_size > min_sample_size || throw(ArgumentError("Excpected max_sample_size > min_sample_size."))
step > 0 || throw(ArgumentError("Excpected step > 0."))
eval_sample_size > max_sample_size || throw(ArgumentError("Excpected eval_sample_size > max_sample_size."))
repeats > 0 || throw(ArgumentError("Excpected repeats > 0."))

# ------------------------------------------------- #

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

sample_sizes = min_sample_size:step:max_sample_size

println("*************************************************")

println("Instance : $(split(instance, "/")[end])")

wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))
pool = kep_pool(wmd_file, dat_file)

perc_cost = zeros(length(sample_sizes), repeats)
real_cost = zeros(length(sample_sizes), repeats)
solve_time = zeros(length(sample_sizes), repeats)

for (s, sample_size) in enumerate(sample_sizes)
    sample_size = Int(sample_size)
    println("  Sample size : $sample_size ...")
    for repeat in 1:repeats
        scenarios = sample_scenarios(pool, generation, sample_size)
        eval_scenarios = sample_scenarios(pool, generation, eval_sample_size)
        v_perc, m_perc = SP(formulation, pool, scenarios, budget)
        v_real, m_real = decision_evaluation(formulation, v_perc["x"], pool, eval_scenarios)
        perc_cost[s,repeat] = m_perc["objective_value"]
        real_cost[s,repeat] = m_real["objective_value"]
        solve_time[s,repeat] = m_perc["solve_time"]
    end
end

mean_perc_cost = mean(perc_cost, dims=2)
var_perc_cost = var(perc_cost, dims=2)
mean_real_cost = mean(real_cost, dims=2)
var_real_cost = var(real_cost, dims=2)
mean_solve_time = mean(solve_time, dims=2)
var_solve_time = var(solve_time, dims=2)

p1 = plot(
    sample_sizes,
    [mean_perc_cost, mean_real_cost],
    label = ["perc cost" "real cost"],
    marker = :circle,
    legend = :bottomright,
    xlabel = "Scenario factor",
    ylabel = "Cost",
)

p2 = plot(
    sample_sizes,
    [var_perc_cost, var_real_cost],
    label = ["perc cost" "real cost"],
    marker = :circle,
    legend = :topright,
    title = split(instance, "/")[end],
    xlabel = "Scenario factor",
    ylabel = "Variance",
)

p3 = plot(
    sample_sizes,
    mean_solve_time,
    label = "",
    marker = :circle,
    legend = :bottomright,
    xlabel = "Scenario factor",
    ylabel = "Solving time (seconds)",
)

plot(p1, p2, p3, layout=(1,3), size=(900,250), bottom_margin=5Plots.mm, left_margin=3Plots.mm)
savefig(joinpath(@__DIR__, joinpath("saves", "scenarios.png")))

println("*************************************************")