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
    "min_budget"
        help = "minimum HLA-crossmatch test budget"
        arg_type = Int
        required = true
    "max_budget"
        help = "maximum HLA-crossmatch test budget"
        arg_type = Int
        required = true
    "step"
        help = "step in HLA-crossmatch budget"
        arg_type = Int
        required = true
    "sample_size"
        help = "number of SAA samples in the SAA"
        arg_type = Int
        required = true
    "eval_sample_size"
        help = "number of SAA samples in the evaluation of a decision"
        arg_type = Int
        required = true
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
min_budget = parsed_args["min_budget"]
max_budget = parsed_args["max_budget"]
step = parsed_args["step"]
sample_size = parsed_args["sample_size"]
eval_sample_size = parsed_args["eval_sample_size"]
maxtime = parsed_args["maxtime"]

min_budget > 0 || throw(ArgumentError("Excpected min_budget > 0."))
max_budget > min_budget || throw(ArgumentError("Excpected max_budget > min_budget."))
step > 0 || throw(ArgumentError("Excpected step > 0."))
sample_size > 0 || throw(ArgumentError("Excpected sample_size > 0."))
eval_sample_size > sample_size || throw(ArgumentError("Excpected eval_sample_size > sample_size."))

# ------------------------------------------------- #

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

budgets = min_budget:step:max_budget

perc_cost = zeros(length(budgets))
real_cost = zeros(length(budgets))
eev_real = zeros(length(budgets))
ws_real = zeros(length(budgets))
solvetime = zeros(length(budgets))

println("*************************************************")
println("Instance : $instance")

wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))
pool = kep_pool(wmd_file, dat_file)

scenarios = sample_scenarios(pool, generation, sample_size)
eval_scenarios = sample_scenarios(pool, generation, eval_sample_size)

for (b, budget) in enumerate(budgets)

    println("  Budget : $budget ...")
    budget = Int(budget)

    v_perc, m_perc = SP(formulation, pool, scenarios, budget, maxtime=maxtime)
    v_real, m_real = decision_evaluation(formulation, v_perc["x"], pool, eval_scenarios, maxtime=maxtime)
    v_evp, m_evp = EVP(formulation, pool, scenarios, budget, maxtime=maxtime)
    v_eev, m_eev = decision_evaluation(formulation, v_evp["x"], pool, eval_scenarios, maxtime=maxtime)
    ws, m_ws = WS(formulation, pool, eval_scenarios, budget, maxtime=maxtime)

    perc_cost[b] = m_perc["objective_value"]
    real_cost[b] = m_real["objective_value"]
    eev_real[b] = m_eev["objective_value"]
    ws_real[b] = ws
    solvetime[b] = m_perc["solve_time"]

end

p1 = plot(
    budgets,
    [perc_cost, real_cost, eev_real, ws_real],
    label = ["perc cost" "real cost" "EEV" "WS"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    title = split(instance, "/")[end],
    xlabel = "Budget factor",
    ylabel = "Objective",
    size = (450, 400)
)

p2 = plot(
    budgets,
    solvetime,
    label = "solve time",
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :topleft,
    title = split(instance, "/")[end],
    xlabel = "Budget factor",
    ylabel = "Seconds",
    size = (450, 400)
)

plot(p1, p2, layout=(1,2), size=(600,250), bottom_margin=5Plots.mm, left_margin=3Plots.mm)
savefig(joinpath(@__DIR__, joinpath("saves", "budget.png")))

println("*************************************************")