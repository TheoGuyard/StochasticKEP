"""
    matching_deterministic
    
Deterministic matching formulation for a KEP problem. Only cycles of length 2 
are allowed by this model.
    
# Parameters
* `pool::MetaDiGraph` : Kep pool
* `scenario::Scenario` : The scenario to consider
* `B::Real` : Maximum number of HLA-crossmatch tests allowed
* `maxtime::Real=60` : Maximum solving time in seconds
* `relax::Bool=false` : Whether to real recourse integrity constraints
* `verbose::Bool=false` : Whether to toogle verbosity
"""
function matching_deterministic(pool::MetaDiGraph, scenario::Scenario, B::Real;
    maxtime::Real=60, relax::Bool=false, verbose::Bool=false)

    assert_pool_scenario_compatibility(pool, scenario)

    matching_pool = extract_matching_pool(pool)
    matching_scenario = extract_matching_scenario(scenario, pool, matching_pool)

    ne(matching_pool) > 0 || throw(ArgumentError("Matching graph obtained from the pool has no edges."))

    V = nv(matching_pool)
    E = ne(matching_pool)
    weights = [get_prop(matching_pool, edge, :weight) for edge in edges(matching_pool)]

    model = Model(optimizer_with_attributes(GLPK.Optimizer, 
        "tm_lim" => round(1000*maxtime), "msg_lev" => GLPK.GLP_MSG_OFF))

    # Variables : 
    #   y_e = 1 if the transplant e is carried out
    #   y_e = 0 else
    if relax
        @variable(model, y[1:E]>=0, base_name="y")
    else
        @variable(model, y[1:E], Bin, base_name="y")
    end

    # Objective : 
    #   Maximize the total transplant utility
    @objective(model, Max, y' * weights)

    # Matching constraint : 
    #   Each vertex is connected to at most one other vertex
    @constraint(model, matching[v=1:V], 
        sum(y[e] for e in num_edge_with_extremity(matching_pool, v)) <= 1)

    # HLA-crossmatch constraint :
    #    Transplantations can only be carried out if the HLA-crossmatch test
    #    permits it  
    @constraint(model, outcome[e=1:E], y[e] <= matching_scenario.outcome[e])

    # Maximum number of HLA-crossmatch tests constraint :
    #   At most B HLA-crossmatch tests can be taken
    @constraint(model, max_number_test, sum(y) <= B/2)

    optimize!(model)

    # Reconstruct solutions for the initial graph
    if termination_status(model) == MOI.OPTIMAL
        y_opt = value.(model[:y])
        y_opt_reconst = reconstruct_matching_variable(y_opt, pool, matching_pool)
    else
        y_opt_reconst = zeros(ne(pool))
    end

    vars = Dict("y" => y_opt_reconst)
    monitoring = get_monitoring(model)

    return vars, monitoring

end

"""
    matching_stochastic

Stochastic matching formulation for a KEP problem. Only cycles of length 2 are 
allowed by this model.
    
# Parameters
* `pool::MetaDiGraph` : Kep pool
* `scenarios::Array{Scenario,1}` : The sample of scenarios
* `B::Real` : Maximum number of HLA-crossmatch tests allowed
* `maxtime::Real=60` : Maximum solving time in seconds
* `relax::Bool=false` : Whether to real recourse integrity constraints
* `verbose::Bool=false` : Whether to toogle verbosity
"""
function matching_stochastic(pool::MetaDiGraph, scenarios::Array{Scenario,1}, 
    B::Real; maxtime::Real=60, relax::Bool=false, verbose::Bool=false,
    fix_x::Union{Nothing,Array{Float64,1}}=nothing)
    
    assert_pool_scenarios_compatibility(pool, scenarios)

    matching_pool = extract_matching_pool(pool)
    matching_scenarios = extract_matching_scenarios(scenarios, pool, matching_pool)

    ne(matching_pool) > 0 || throw(ArgumentError("Matching graph obtained from the pool has no edges."))

    model = Model(optimizer_with_attributes(GLPK.Optimizer, 
        "tm_lim" => round(1000*maxtime), "msg_lev" => GLPK.GLP_MSG_OFF))

    K = length(scenarios)
    V = nv(matching_pool)
    E = ne(matching_pool)
    weights = [get_prop(matching_pool, edge, :weight) for edge in edges(matching_pool)]
    probas = [scenario.proba for scenario in matching_scenarios]

    # First-stage variables : 
    #   x_e = 1 if the transplant e takes a HLA-crossmatch test
    #   x_e = 0 else
    @variable(model, x[1:E], Bin, base_name="x")

    # Second-stage variables : 
    #   y_e^k = 1 if the transplant e is carried out in the scenario k
    #   y_e^k = 0 else
    if relax
        @variable(model, 1>=y[1:E,1:K]>=0, base_name="y")
    else
        @variable(model, y[1:E,1:K], Bin, base_name="y")
    end

    # Objective : 
    #   Maximize the total transplant utility over the scenarios
    @objective(model, Max, (y' * weights)' * probas)

    # Matching constraint : 
    #   Each vertex is connected to at most one other vertex
    @constraint(model, matching[v=1:V, k=1:K], 
        sum(y[e,k] for e in num_edge_with_extremity(matching_pool, v)) <= 1)

    # Testing and outcome constraint :
    #   Only tested transplants with a negative HLA-crosslatch outcomes can be 
    #   carried out
    @constraint(model, tested[e=1:E, k=1:K], y[e,k] <= x[e])
    @constraint(model, outcome[e=1:E, k=1:K], 
        y[e,k] <= matching_scenarios[k].outcome[e])

    # Maximum number of HLA-crossmatch tests constraint :
    #   At most B HLA-crossmatch tests can be taken
    @constraint(model, max_number_test, sum(x) <= B/2)

    # Fix decision variables if needed
    if !(fix_x isa Nothing)
        fix_x = get_fix_matching(fix_x, pool, matching_pool)
        length(model[:x]) == length(fix_x) || throw("Dimension of fix_x and of the variable x missmatch.")
        for e in 1:E
            fix(model[:x][e], fix_x[e]; force = true)
        end
    end

    optimize!(model)

    # Reconstruct solutions for the initial graph
    if termination_status(model) == MOI.OPTIMAL
        x_opt = value.(model[:x])
        y_opt = value.(model[:y])
        x_opt_reconst = reconstruct_matching_variable(x_opt, pool, matching_pool)
        y_opt_reconst = zeros(ne(pool),K)
        for k in 1:K
            y_opt_reconst[:,k] = reconstruct_matching_variable(y_opt[:,k], pool, matching_pool)
        end
    else
        x_opt_reconst = zeros(ne(pool))
        y_opt_reconst = zeros(ne(pool),K)
    end

    vars = Dict("x" => x_opt_reconst, "y" => y_opt_reconst)
    monitoring = get_monitoring(model)

    return vars, monitoring
    
end