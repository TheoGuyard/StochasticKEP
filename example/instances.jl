using DataFrames
using LightGraphs
using MetaGraphs

include("../src/StochasticKep.jl")
using .StochasticKep

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")

instances = [
    "preflib-md/MD-00001-00000075",
    "preflib-md/MD-00001-00000080",
    "preflib-md/MD-00001-00000085",
    "preflib-md/MD-00001-00000090",
    "preflib-md/MD-00001-00000095",
    "preflib-md/MD-00001-00000100",
    "preflib-md/MD-00001-00000105",
]

budget_factor = 1
scenario_factor = 1
formulation = "matching"
generation = "BinomialUNOS"

table = DataFrame(
    instance = String[], 
    formulation = String[],
    V = Int[],
    Va = Int[],
    A = Int[],
    generation = String[],
    budget = Real[],
    nbscenarios = Int[],
    card_x = Int[],
    card_y = Int[]
)

for instance in instances

    wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
    dat_file = joinpath(data_folder, join([instance, ".dat"]))
    pool = kep_pool(wmd_file, dat_file)

    if formulation == "matching"
        matching_pool = extract_matching_pool(pool)
        removed_vertices = []
        removed_a_vertices = []
        for v in 1:nv(matching_pool)
            if length(neighbors(matching_pool, v)) == 0
                push!(removed_vertices, v)
                if get_prop(matching_pool, v, :altruist)
                    push!(removed_a_vertices, v)
                end
            end
        end
        vertices = nv(matching_pool) - length(removed_vertices)
        a_vertices = get_prop(matching_pool, :nb_altruist) - length(removed_a_vertices)
        arcs = ne(matching_pool)
    else
        vertices = nv(pool)
        a_vertices = get_prop(pool, :nb_altruist)
        arcs = ne(pool) 
    end

    budget = convert(Int, ceil(budget_factor * vertices))
    nb_scenarios = convert(Int, ceil(scenario_factor * arcs))

    push!(table, [
        instance,
        formulation,
        vertices,
        a_vertices,
        arcs,
        generation,
        budget,
        nb_scenarios,
        arcs,
        arcs * nb_scenarios,
    ])
end

print(table)