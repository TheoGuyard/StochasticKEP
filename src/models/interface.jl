function DP(formulation::String, pool::MetaDiGraph, scenario::Scenario, B::Real;
    maxtime::Real=60, relax::Bool=false, verbose::Bool=false)

    formulation in FORMULATIONS || throw(ArgumentError("Formulation not understood. Allowed formulations : $FORMULATIONS."))
    
    if formulation == "matching"
        return matching_deterministic(pool, scenario, B, maxtime=maxtime, 
            relax=relax, verbose=verbose)
    elseif formulation == "relaxed_arc"
        return relaxed_arc_deterministic(pool, scenario, B, maxtime=maxtime, 
        relax=relax, verbose=verbose)
    end

end

function SP(formulation::String, pool::MetaDiGraph, scenarios::Array{Scenario,1}, 
    B::Real; maxtime::Real=60, relax::Bool=false, verbose::Bool=false, 
    fix_x::Union{Nothing,Array{Float64,1}}=nothing)
    
    formulation in FORMULATIONS || throw(ArgumentError("Formulation not understood. Allowed formulations : $FORMULATIONS."))
    
    if formulation == "matching"
        return matching_stochastic(pool, scenarios, B, maxtime=maxtime, 
        relax=relax, verbose=verbose, fix_x=fix_x)
    elseif formulation == "relaxed_arc"
        return relaxed_arc_stochastic(pool, scenarios, B, maxtime=maxtime, 
        relax=relax, verbose=verbose, fix_x=fix_x)
    end
end

function WS(formulation::String, pool::MetaDiGraph, scenarios::Array{Scenario,1},
    B::Real; maxtime::Real=60, verbose::Bool=false)

    formulation in FORMULATIONS || throw(ArgumentError("Formulation not understood. Allowed formulations : $FORMULATIONS."))

    time_left = maxtime

    ws = 0.
    for (s, scenario) in enumerate(scenarios)
        start_time = time()
        vars, monitoring = DP(formulation, pool, scenario, B, maxtime=time_left, verbose=verbose)
        if monitoring["termination_status"] == MOI.OPTIMAL
            ws += monitoring["objective_value"]
        else
            throw(ErrorException("Wait-and-see problem for scenario $s not solved."))
            ws -= Inf
        end
        time_left -= time() - start_time
    end

    ws /= length(scenarios)

    monitoring = Dict(
        "termination_status" => ws != -Inf ? MOI.OPTIMAL : MOI.TIME_LIMIT,
        "solve_time" => (maxtime - time_left)
    )

    return ws, monitoring

end

function EVP(formulation::String, pool::MetaDiGraph, scenarios::Array{Scenario,1},
    B::Real; maxtime::Real=60, verbose::Bool=false)

    formulation in FORMULATIONS || throw(ArgumentError("Formulation not understood. Allowed formulations : $FORMULATIONS."))

    return SP(formulation, pool, [expected_scenario(scenarios)], B,
        maxtime=maxtime, relax=true, verbose=verbose)

end

function decision_evaluation(formulation::String, x::Array{Float64,1}, 
    pool::MetaDiGraph, scenarios::Array{Scenario}; maxtime::Real=60, 
    relax::Bool=false, verbose::Bool=false)

    formulation in FORMULATIONS || throw(ArgumentError("Formulation not understood. Allowed formulations : $FORMULATIONS."))
    maximum(x-round.(x)) < 1e-4 || relax || throw(ArgumentError("Decision x is not integer. Need to set relax=true."))

    dummy_B = convert(Int, ceil(sum(x)))
    vars, monitoring = SP(formulation, pool, scenarios, dummy_B, 
        maxtime=maxtime, relax=relax, verbose=verbose, fix_x=x)
    
    return vars, monitoring

end

function find_cycles(y::Array{Float64,1}, pool::MetaDiGraph; maxlength::Int=4)

    all(x -> max(abs.(x-round(x))) < 1e-4, y) || throw(ArgumentError("y has some fractionnal coefficents."))
    length(y) == ne(pool) || throw(ArgumentError("Size of y missmatch with number of arcs in the pool."))

    reduced_pool = copy(pool)
    for (a, arc) in enumerate(edges(pool))
        if abs(y[a]) < 1e-4
            rem_edge!(reduced_pool, arc)
        end
    end
   
    cycles = simplecycles(reduced_pool)

    nb_cycles = Dict()
    for i in 1:(maxlength-1)
        nb_cycles["$i"] = length([cycle for cycle in cycles if length(cycle) == i])
    end
    nb_cycles["$(maxlength)+"] = length([cycle for cycle in cycles if length(cycle) >= maxlength])

    return cycles, nb_cycles
end

function find_cycles(y::Array{Float64,2}, pool::MetaDiGraph; maxlength::Int=4)

    K = size(y)[2]

    cycles = []
    mean_nb_cycles = Dict()
    for i in 1:(maxlength-1)
        mean_nb_cycles["$i"] = 0
    end
    mean_nb_cycles["$(maxlength)+"] = 0

    for k in 1:K
        scenario_cycles, scenario_nb_cycles = find_cycles(y[:,k], pool, maxlength=maxlength)
        push!(cycles, scenario_cycles)
        for key in keys(mean_nb_cycles)
            mean_nb_cycles[key] += scenario_nb_cycles[key]
        end
    end

    for key in keys(mean_nb_cycles)
        mean_nb_cycles[key] /= K
    end

    return cycles, mean_nb_cycles

end