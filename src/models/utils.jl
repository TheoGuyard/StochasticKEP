"""
    extract_matching_pool

Return an undirected graph with the same vertex as the one in the original 
pool. An edge i-j is created iff there exsits a directed edges i->j and j->i in 
orinial pool. The new weight is the sum of the two weights in orinial pool.

# Parameters
* `pool::MetaDiGraph` : Orinigal KEP pool
"""
function extract_matching_pool(pool::MetaDiGraph)

    # The extracted graph is no longer a directed graph
    matching_pool = MetaGraph(nv(pool), 0)

    # Copy metadata
    set_prop!(matching_pool, :nb_altruist, get_prop(pool, :nb_altruist))
    for i in vertices(matching_pool)
        set_props!(matching_pool, i, props(pool, i))
    end

    # Create edges i-j iif i->j and j->i exists in pool and set the metadata
    for i in vertices(matching_pool), j in vertices(matching_pool)
        if has_edge(pool, i, j) & has_edge(pool, j, i)
            add_edge!(matching_pool, i, j)
            set_props!(matching_pool, i, j, Dict(
                :weight => get_prop(pool, i, j, :weight) +
                           get_prop(pool, j, i, :weight),
            ))
        end
    end

    return matching_pool

end

"""
    extract_matching_scenario

Return the merged scenario of a matching formulation. The weigth of i-j in the 
matching formulation is the product of the weights of i->j and j->i in the 
original pool.

# Parameters
* `scenario::Scenario` : Original scenario
* `pool::MetaDiGraph` : Original KEP pool
* `matching_pool::MetaGraph` : Matching KEP pool
"""
function extract_matching_scenario(scenario::Scenario, pool::MetaDiGraph, matching_pool::MetaGraph)
    matching_outcome = zeros(ne(matching_pool))
    for e in 1:ne(matching_pool)
        matching_edge = collect(edges(matching_pool))[e]
        e1 = num_edge_with_extremity(pool, matching_edge.src, matching_edge.dst)
        e2 = num_edge_with_extremity(pool, matching_edge.dst, matching_edge.src)
        matching_outcome[e] = scenario.outcome[e1] * scenario.outcome[e2]
    end
    return Scenario(scenario.proba, scenario.generation, scenario.failure_rate, matching_outcome)
end

"""
    extract_matching_scenarios

Same as `extract_matching_scenario` but with a sample of scenarios.

# Parameters
* `scenarios::Array{Scenario,1}` : Original scenario sample
* `pool::MetaDiGraph` : Original KEP pool
* `matching_pool::MetaGraph` : Matching KEP pool
"""
function extract_matching_scenarios(scenarios::Array{Scenario,1}, 
    pool::MetaDiGraph, matching_pool::MetaGraph)

    matching_scenarios = []

    for scenario in scenarios
        push!(matching_scenarios, extract_matching_scenario(scenario, pool, matching_pool))
    end

    return matching_scenarios

end

function reconstruct_matching_variable(var::Array{Float64,1}, 
    pool::MetaDiGraph, matching_pool::MetaGraph)

    var_reconst = zeros(ne(pool))

    for (e, edge) in enumerate(edges(matching_pool))
        source = edge.src
        destination = edge.dst
        has_edge(pool, source, destination) || throw(ArgumentError("Edge $source => $destination found in the matching pool is not in the original pool."))
        has_edge(pool, destination, source) || throw(ArgumentError("Edge $destination => $source found in the matching pool is not in the original pool."))
        if var[e] != 0
            arc_1 = num_edge_with_extremity(pool, source, destination)
            arc_2 = num_edge_with_extremity(pool, destination, source)
            var_reconst[arc_1] = var[e]
            var_reconst[arc_2] = var[e]
        end
    end

    return var_reconst
end

function get_fix_matching(fix::Array{Float64,1}, pool::MetaDiGraph, 
    matching_pool::MetaGraph)

    fix_matching = zeros(ne(matching_pool))

    for (e, edge) in enumerate(edges(matching_pool))
        source = edge.src
        destination = edge.dst
        has_edge(pool, source, destination) || throw(ArgumentError("Edge $source => $destination found in the matching pool is not in the original pool."))
        has_edge(pool, destination, source) || throw(ArgumentError("Edge $destination => $source found in the matching pool is not in the original pool."))
        fix_1 = fix[num_edge_with_extremity(pool, source, destination)]
        fix_2 = fix[num_edge_with_extremity(pool, destination, source)]
        fix_1 == fix_2 || throw(ErrorException("Different fixed values for the same matching edge."))
        fix_matching[e] = fix_1
    end

    return fix_matching
    
end

function get_monitoring(model::JuMP.Model)
    monitoring = Dict(
        "termination_status" => try termination_status(model) catch; MOI.OPTIMIZE_NOT_CALLED end,
        "objective_value" => try objective_value(model) catch; -Inf end,
        "solve_time" => try solve_time(model) catch; -1 end,
        "objective_bound" => try objective_bound(model) catch; Inf end,
        "relative_gap" => try relative_gap(model) catch; Inf end,
        "simplex_iterations" => try simplex_iterations(model) catch; -1 end,
        "barrier_iterations" => try barrier_iterations(model) catch; -1 end,
        "node_count" => try node_count(model) catch; -1 end,
    )
    return monitoring
end

function matching_budget(matching_pool::MetaGraph, budget_factor::Real)

    not_connected_vertex = 0
    for v in 1:nv(matching_pool)
        if length(neighbors(matching_pool, v)) == 0
            not_connected_vertex += 1
        end
    end

    budget = budget_factor * (nv(matching_pool) - not_connected_vertex)

    return budget

end

function num_edge_with_extremity(graph::Union{MetaDiGraph,MetaGraph}, v::Int)
    has_vertex(graph, v) || throw(ArgumentError("Vertex $v is not in the graph"))
    res = [e for (e, edge) in enumerate(edges(graph)) if 
        (edge.src == v) | (edge.dst == v)]
    return res
end

function num_edge_with_extremity(graph::Union{MetaDiGraph,MetaGraph}, vs::Int, vd::Int)
    has_vertex(graph, vs) || error("Vertex $vs is not in the graph")
    has_vertex(graph, vd) || error("Vertex $vd is not in the graph")
    if graph isa MetaDiGraph
        res = [e for (e, edge) in enumerate(edges(graph)) if 
            (edge.src == vs) & (edge.dst == vd)]
    elseif graph isa MetaGraph
        res = [e for (e, edge) in enumerate(edges(graph)) if 
            ((edge.src == vs) & (edge.dst == vd)) | 
            ((edge.src == vd) & (edge.dst == vs))]
    end
    return res[1]
end

function num_edge_with_source(graph::MetaDiGraph, v::Int)
    has_vertex(graph, v) || throw(ArgumentError("Vertex $v is not in the graph"))
    res = [e for (e, edge) in enumerate(edges(graph)) if (edge.src == v)]
    return res
end

function num_edge_with_destination(graph::MetaDiGraph, v::Int)
    has_vertex(graph, v) || throw(ArgumentError("Vertex $v is not in the graph"))
    res = [e for (e, edge) in enumerate(edges(graph)) if (edge.dst == v)]
    return res
end

function assert_pool_scenario_compatibility(pool::MetaDiGraph, scenario::Scenario)
    ne(pool) == length(scenario.outcome) || throw(ArgumentError("Length of scenario outcomes doesn't match the number of transplants in the KEP pool."))
end

function assert_pool_scenarios_compatibility(pool::MetaDiGraph, scenarios::Array{Scenario,1})
    length(scenarios) > 0 || throw(ArgumentError("Empty scenario array."))
    for scenario in scenarios
        assert_pool_scenario_compatibility(pool, scenario)
    end
end

