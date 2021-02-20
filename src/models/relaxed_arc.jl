"""
    relaxed_arc_deterministic
    
Deterministic relaxed-arc formulation for a KEP problem. No cycle length 
limitation is imposed.
    
# Parameters
* `pool::MetaDiGraph` : Kep pool
* `scenario::Scenario` : The scenario to consider
* `B::Real` : Maximum number of HLA-crossmatch tests allowed
* `maxtime::Real=60` : Maximum solving time in seconds
* `relax::Bool=false` : Whether to relax recourse integrity constraints
* `verbose::Bool=false` : Whether to toogle verbosity
"""
function relaxed_arc_deterministic(pool::MetaDiGraph, scenario::Scenario, B::Real; 
    maxtime::Real=60, relax::Bool=false, verbose::Bool=false)

    assert_pool_scenario_compatibility(pool, scenario)

    V = nv(pool)
    A = ne(pool)
    weights = [get_prop(pool, arc, :weight) for arc in edges(pool)]

    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
        "TimeLimit" => maxtime, "OutputFlag" => 0))

    # Variables : 
    #   x_a = 1 if the transplant a is carried out
    #   x_a = 0 else
    if relax
        @variable(model, x[1:A] >= 0, base_name="x")
    else
        @variable(model, x[1:A], Bin, base_name="x")
    end

    # Objective : 
    #   Maximize the total transplant utility
    @objective(model, Max, x' * weights)
    
    # Flow constraint 1 : 
    #   Each vertex has more entering edges than leaving edges
    @constraint(model, flow1[v=1:V], 
        sum(x[a] for a in num_edge_with_destination(pool, v)) <= 
        sum(x[a] for a in num_edge_with_source(pool, v))
    )
    
    # Flow constraint 2 : 
    #   Each vertex has at most one leaving edge
    @constraint(model, flow2[v=1:V], 
        sum(x[a] for a in num_edge_with_destination(pool, v)) <= 1)

    # Outcome constraint :
    #    Transplantations can only be carried out if the HLA-crossmatch test
    #    permits it  
    @constraint(model, outcome[a=1:A], x[a] <= scenario.outcome[a])

    # Maximum number of HLA-crossmatch tests constraint :
    #   In total, B HLA-crossmatch tests can be taken to test whether a
    #   transplant can be carried out or not
    @constraint(model, max_number_test, sum(x) <= B)

    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        x_opt = value.(model[:x])
    else
        x_opt = zeros(A)
    end

    vars = Dict("x" => x_opt)
    monitoring = get_monitoring(model)

    return vars, monitoring

end

"""
    relaxed_arc_stochastic

Stochastic relaxed-arc formulation for a KEP problem. No cycle length 
limitation is imposed.
    
# Parameters
* `pool::MetaDiGraph` : Kep pool
* `scenarios::Array{Scenario,1}` : The sample of scenarios
* `B::Real` : Maximum number of HLA-crossmatch tests allowed
* `maxtime::Real=60` : Maximum solving time in seconds
* `relax::Bool=false` : Whether to relax recourse integrity constraints
* `verbose::Bool=false` : Whether to toogle verbosity
"""
function relaxed_arc_stochastic(pool::MetaDiGraph, scenarios::Array{Scenario,1}, 
    B::Real; maxtime::Real=60, relax::Bool=false, verbose::Bool=false,
    fix_x::Union{Nothing,Array{Float64,1}}=nothing)

    assert_pool_scenarios_compatibility(pool, scenarios)

    model = Model(optimizer_with_attributes(Gurobi.Optimizer, 
        "TimeLimit" => maxtime, "OutputFlag" => 0))

    K = length(scenarios)
    V = nv(pool)
    A = ne(pool)
    weights = [get_prop(pool, arc, :weight) for arc in edges(pool)]
    probas = [scenario.proba for scenario in scenarios]
    
    # First-stage variables : 
    #   x_a = 1 if the transplant a takes a HLA-crossmatch test
    #   x_a = 0 else
    @variable(model, x[1:A], Bin, base_name="x")

    # Second-stage variables : 
    #   y_a^k = 1 if the transplant e is carried out in the scenario k
    #   y_a^k = 0
    if relax
        @variable(model, y[1:A,1:K]>=0, base_name="y")
    else
        @variable(model, y[1:A,1:K], Bin, base_name="y")
    end

    # Objective : 
    #   Maximize the total transplant utility over the scenarios
    @objective(model, Max, (y' * weights)' * probas)

    # Flow constraint 1 : 
    #   Each vertex has more entering edges than leaving edges
    @constraint(model, flow1[v=1:V, k=1:K], 
        sum(y[a,k] for a in num_edge_with_destination(pool, v)) <= 
        sum(y[a,k] for a in num_edge_with_source(pool, v))
        )

    # Flow constraint 2 : 
    #   Each vertex has at most one leaving edge
    @constraint(model, flow2[v=1:V, k=1:K], 
        sum(y[a,k] for a in num_edge_with_destination(pool, v)) <= 1)

    # Maximum number of HLA-crossmatch tests constraint :
    #   In total, B HLA-crossmatch tests can be taken to test whether a
    #   transplant can be carried out or not
    @constraint(model, max_number_test, sum(x) <= B)

    # Testing and outcome constraint :
    #   Only tested transplants with a negative HLA-crosslatch outcomes can be 
    #   carried out
    @constraint(model, tested[a=1:A, k=1:K], y[a,k] <= x[a])
    @constraint(model, outcome[a=1:A, k=1:K], 
        y[a,k] <= scenarios[k].outcome[a])

    # Fix decision variables if needed
    if !(fix_x isa Nothing)
        length(model[:x]) == length(fix_x) || throw("Dimension of fix_x and of the variable x missmatch.")
        for a in 1:A
            fix(model[:x][a], fix_x[a]; force = true)
        end
    end

    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        x_opt = value.(model[:x])
        y_opt = value.(model[:y])
    else
        x_opt = zeros(A)
        y_opt = zeros(A,K)
    end

    vars = Dict("x" => x_opt, "y" => y_opt)
    monitoring = get_monitoring(model)

    return vars, monitoring

end