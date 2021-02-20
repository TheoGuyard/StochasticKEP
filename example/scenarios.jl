using DataFrames
using LightGraphs
using Plots
using Statistics

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

instance = "preflib-md/MD-00001-00000040"

scenario_factors = 0.2:0.1:1.2
repeats = 30

formulation = "matching"
generation = "BinomialUNOS"
budget_factor = 1

println("\nInstance : $(split(instance, "/")[end])")

wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))
pool = kep_pool(wmd_file, dat_file)

perc_cost = zeros(length(scenario_factors), repeats)
real_cost = zeros(length(scenario_factors), repeats)
solve_time = zeros(length(scenario_factors), repeats)

for (f, factor) in enumerate(scenario_factors)
    println("  Factor : $factor")

    if formulation == "matching"
        matching_pool = extract_matching_pool(pool)
        nb_scenarios = convert(Int, ceil(factor * ne(matching_pool)))
        budget = convert(Int, ceil(matching_budget(matching_pool, budget_factor)))
    else
        nb_scenarios = convert(Int, ceil(factor * ne(pool)))
        budget = convert(Int, ceil(budget_factor * nv(pool)))
    end

    for repeat in 1:repeats
        scenarios = sample_scenarios(pool, generation, nb_scenarios)
        eval_scenarios = sample_scenarios(pool, generation, 5000)
        v_perc, m_perc = SP(formulation, pool, scenarios, budget)
        v_real, m_real = decision_evaluation(formulation, v_perc["x"], pool, eval_scenarios)
        perc_cost[f,repeat] = m_perc["objective_value"]
        real_cost[f,repeat] = m_real["objective_value"]
        solve_time[f,repeat] = m_perc["solve_time"]
    end
end

mean_perc_cost = mean(perc_cost, dims=2)
var_perc_cost = var(perc_cost, dims=2)
mean_real_cost = mean(real_cost, dims=2)
var_real_cost = var(real_cost, dims=2)
mean_solve_time = mean(solve_time, dims=2)
var_solve_time = var(solve_time, dims=2)

p1 = plot(
    scenario_factors,
    [mean_perc_cost, mean_real_cost],
    label = ["perc cost" "real cost"],
    marker = :circle,
    legend = :bottomright,
    xlabel = "Scenario factor",
    ylabel = "Cost",
)

p2 = plot(
    scenario_factors,
    [var_perc_cost, var_real_cost],
    label = ["perc cost" "real cost"],
    marker = :circle,
    legend = :topright,
    title = split(instance, "/")[end],
    xlabel = "Scenario factor",
    ylabel = "Variance",
)

p3 = plot(
    scenario_factors,
    mean_solve_time,
    label = "",
    marker = :circle,
    legend = :bottomright,
    xlabel = "Scenario factor",
    ylabel = "Solving time (seconds)",
)

plot(p1, p2, p3, layout=(1,3), size=(900,250), bottom_margin=5Plots.mm, left_margin=3Plots.mm)