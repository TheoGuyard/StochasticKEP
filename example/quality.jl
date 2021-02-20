using DataFrames
using LightGraphs
using Plots

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

instances = [
    "preflib-md/MD-00001-00000001",
    "preflib-md/MD-00001-00000003",
    "preflib-md/MD-00001-00000005",
    "preflib-md/MD-00001-00000007",
    "preflib-md/MD-00001-00000009",
    "preflib-md/MD-00001-00000011",
    "preflib-md/MD-00001-00000013",
    "preflib-md/MD-00001-00000015",
    "preflib-md/MD-00001-00000017",
    "preflib-md/MD-00001-00000019",
]

sp_perc = []
sp_real = []
spr_perc = []
spr_real = []
ev_perc = []
ev_real = []
eev_perc = []
eev_real = []
ws_values = []

formulation = "matching"
generation = "BinomialUNOS"

scenario_factor = 10
budget_factor = (2/3)
evaluation_factor = 10

for instance in instances

    wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
    dat_file = joinpath(data_folder, join([instance, ".dat"]))
    pool = kep_pool(wmd_file, dat_file)

    if formulation == "matching"
        matching_pool = extract_matching_pool(pool)
        nb_scenarios = convert(Int, ceil(scenario_factor * ne(matching_pool)))
        nb_scenarios_eval = convert(Int, ceil(evaluation_factor * ne(matching_pool)))
        budget = convert(Int, ceil(budget_factor * ne(matching_pool)))
    else
        nb_scenarios = convert(Int, ceil(scenario_factor * ne(pool)))
        nb_scenarios_eval = convert(Int, ceil(evaluation_factor * ne(pool)))
        budget = convert(Int, ceil(budget_factor * ne(pool)))
    end

    println("\nInstance\t: $(split(instance, "/")[end])")
    println("Scenarios\t: $(nb_scenarios)")
    println("Budget\t\t: $(budget)")

    scenarios = sample_scenarios(pool, generation, nb_scenarios)
    eval_scenarios = sample_scenarios(pool, generation, nb_scenarios_eval)

    println("Solving the stochastic problem ...")
    v_sp, m_sp = SP(formulation, pool, scenarios, budget)
    println("  Solve time : $(round(m_sp["solve_time"], digits=2)) seconds")
    push!(sp_perc, m_sp["objective_value"])

    println("Evaluate the decision on 10x more scenarios ...")
    v_sp_eval, m_sp_eval = decision_evaluation(formulation, v_sp["x"], pool, eval_scenarios)
    println("  Evaluation time : $(round(m_sp_eval["solve_time"], digits=2)) seconds")
    push!(sp_real, m_sp_eval["objective_value"])

    println("Solving the stochastic problem (relaxed) ...")
    v_spr, m_spr = SP(formulation, pool, scenarios, budget, relax=true)
    println("  Solve time : $(round(m_spr["solve_time"], digits=2)) seconds")
    push!(spr_perc, m_spr["objective_value"])

    println("Evaluate the decision on 10x more scenarios ...")
    v_spr_eval, m_spr_eval = decision_evaluation(formulation, v_spr["x"], pool, eval_scenarios, relax=true)
    println("  Evaluation time : $(round(m_spr_eval["solve_time"], digits=2)) seconds")
    push!(spr_real, m_spr_eval["objective_value"])

    println("Solving the EEV problem ...")
    v_ev, m_ev, v_eev, m_eev = EEV(formulation, pool, scenarios, budget)
    println("  Solve time : $(round(m_ev["solve_time"] + m_eev["solve_time"], digits=2)) seconds")
    push!(ev_perc, m_ev["objective_value"])
    push!(eev_perc, m_eev["objective_value"])

    println("Evaluate the decision on 10x more scenarios ...")
    v_ev_eval, m_ev_eval = decision_evaluation(formulation, v_ev["x"], pool, eval_scenarios, relax=true)
    push!(ev_real, m_ev_eval["objective_value"])
    println("  Evaluation time : $(round(m_ev_eval["solve_time"], digits=2)) seconds")
    v_eev_eval, m_eev_eval = decision_evaluation(formulation, v_eev["x"], pool, eval_scenarios, relax=true)
    println("  Evaluation time : $(round(m_eev_eval["solve_time"], digits=2)) seconds")
    push!(eev_real, m_eev_eval["objective_value"])

    println("Solving the WS problem ...")
    ws, m_ws = WS(formulation, pool, eval_scenarios)
    println("  Solve time : $(round(m_ws["solve_time"])) seconds")
    push!(ws_values, ws)

end

df = DataFrame(
    sp_perc = Float64[],
    sp_real = Float64[],
    spr_perc = Float64[],
    spr_real = Float64[],
    ev_perc = Float64[],
    ev_real = Float64[],
    eev_perc = Float64[],
    eev_real = Float64[],
    ws_values = Float64[],
    )

for i in 1:length(instances)
    push!(df, [
        sp_perc[i], 
        sp_real[i], 
        spr_perc[i], 
        spr_real[i], 
        ev_perc[i], 
        ev_real[i], 
        eev_perc[i], 
        eev_real[i], 
        ws_values[i]
        ]
    )
end

print(df)

p1 = plot(
    1:length(instances), 
    [100*(sp_real-sp_perc)./sp_real 100*(spr_real-spr_perc)./spr_real 100*(ev_real-ev_perc)./ev_real 100*(eev_real-eev_perc)./eev_real], 
    label = ["sp" "spr" "ev" "eev"], 
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (900, 400)
    )

p2 = plot(
    1:length(instances), 
    [sp_perc-spr_perc sp_perc-eev_perc sp_perc-ws_values sp_perc-ev_perc], 
    label = ["sp_perc-spr_perc" "sp_perc-eev_perc" "sp_perc-ws_values" "sp_perc-ev_perc"], 
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (900, 400)
    )

p3 = plot(
    1:length(instances), 
    [sp_real-spr_real sp_real-eev_real sp_real-ws_values sp_real-ev_real], 
    label = ["sp_real-spr_real" "sp_real-eev_real" "sp_real-ws_values" "sp_real-ev_real"], 
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (900, 400)
    )

plot(p1, p2, p3, layout=(3,1))