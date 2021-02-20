using DataFrames
using LightGraphs
using Plots
using Statistics

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

instance = "preflib-md/MD-00001-00000080"

budget_factors = [0.1,0.2,0.4,0.6,0.8,1,1.2]

perc_cost = zeros(length(budget_factors))
real_cost = zeros(length(budget_factors))
eev_real = zeros(length(budget_factors))
ws_real = zeros(length(budget_factors))
solvetime = zeros(length(budget_factors))

formulation = "matching"
generation = "BinomialUNOS"

scenario_factor = 1
evaluation_factor = 10
maxtime = 600

println("\nInstance : $(split(instance, "/")[end])")

wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
dat_file = joinpath(data_folder, join([instance, ".dat"]))
pool = kep_pool(wmd_file, dat_file)

if formulation == "matching"
    matching_pool = extract_matching_pool(pool)
    nb_scenarios = convert(Int, ceil(scenario_factor * ne(matching_pool)))
    nb_scenarios_eval = convert(Int, ceil(evaluation_factor * ne(matching_pool)))
else
    nb_scenarios = convert(Int, ceil(scenario_factor * ne(pool)))
    nb_scenarios_eval = convert(Int, ceil(evaluation_factor * ne(pool)))
end

scenarios = sample_scenarios(pool, generation, nb_scenarios)
eval_scenarios = sample_scenarios(pool, generation, nb_scenarios_eval)

ref_scenario = Scenario(pool, 1, "NoFailure")
ref_var, ref_monitor = DP(formulation, pool, ref_scenario, ne(pool), maxtime=maxtime)
z_max = ref_monitor["objective_value"]

for (b, budget_factor) in enumerate(budget_factors)

    println("  Budget factor : $(budget_factor)")

    budget = convert(Int, ceil(budget_factor * z_max))

    v_perc, m_perc = SP(formulation, pool, scenarios, budget, maxtime=maxtime)
    v_real, m_real = decision_evaluation(formulation, v_perc["x"], pool, eval_scenarios, maxtime=maxtime)
    v_ev, m_ev, v_eev, m_eev = EEV(formulation, pool, scenarios, budget, maxtime=maxtime)
    ws, m_ws = WS(formulation, pool, eval_scenarios, budget, maxtime=maxtime)

    perc_cost[b] = m_perc["objective_value"]
    real_cost[b] = m_real["objective_value"]
    eev_real[b] = m_eev["objective_value"]
    ws_real[b] = ws
    solvetime[b] = m_perc["solve_time"]

end

p1 = plot(
    budget_factors,
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
    budget_factors,
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