using DataFrames
using LightGraphs
using Plots

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

formulation = "relaxed_arc"
generation = "BinomialUNOS"

scenario_factor = 1
budget_factor = (2/3)

maxlength = 4

colnames = ["$i" for i in 1:(maxlength-1)]
push!(colnames, "$(maxlength)+")
nb_cycles = DataFrame(zeros(length(instances), maxlength), Symbol.(colnames))

for (i, instance) in enumerate(instances)

    wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
    dat_file = joinpath(data_folder, join([instance, ".dat"]))
    pool = kep_pool(wmd_file, dat_file)

    if formulation == "matching"
        matching_pool = extract_matching_pool(pool)
        nb_scenarios = convert(Int, ceil(scenario_factor * ne(matching_pool)))
        budget = convert(Int, ceil(budget_factor * ne(matching_pool)))
    else
        nb_scenarios = convert(Int, ceil(scenario_factor * ne(pool)))
        budget = convert(Int, ceil(budget_factor * ne(pool)))
    end

    println("\nInstance\t: $(split(instance, "/")[end])")
    println("Scenarios\t: $(nb_scenarios)")
    println("Budget\t\t: $(budget)")

    scenarios = sample_scenarios(pool, generation, nb_scenarios)

    println("Solving the stochastic problem ...")
    v_sp, m_sp = SP(formulation, pool, scenarios, budget)
    
    cycles, mean_nb_cycles = find_cycles(v_sp["y"], pool, maxlength=maxlength)
    println(mean_nb_cycles)
    for (k, key) in enumerate(keys(mean_nb_cycles))
        nb_cycles[i,key] = mean_nb_cycles[key]
    end
end

print(nb_cycles)

plot(
    1:length(instances), 
    convert(Matrix,nb_cycles), 
    label = permutedims(names(nb_cycles)),
    seriestype = :scatter, 
    markershape = :circle,
    markersize = 10,
    legend = :bottomright,
    xlabel = "Instance index",
    ylabel = "Mean number of cycles per length",
    size = (900, 400)
    )