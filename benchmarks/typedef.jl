struct BenchmarkInstance
    pool::MetaDiGraph
    formulation::String
    budget::Real
    generation::String
    K::Int 
    K_eval::Int 
    scenarios::Array{Scenario,1}
    eval_scenarios::Array{Scenario,1}
    maxtime::Real
    maxcyclelength::Int
    results::Dict{String,Any}
end

struct BudgetInstance
    pool::MetaDiGraph
    formulation::String
    generation::String
    nb_scenarios::Int
    nb_scenarios_eval::Int
    scenarios::Array{Scenario,1}
    eval_scenarios::Array{Scenario,1}
    maxtime::Real
    budget_factors::Array{Float64,1}
    z_sp_perc::Array{Float64,1}
    z_sp_real::Array{Float64,1}
    z_spr_perc::Array{Float64,1}
    z_spr_real::Array{Float64,1}
    eev::Array{Float64,1}
    ws::Array{Float64,1}
    t_sp_perc::Array{Float64,1}
    t_sp_real::Array{Float64,1}
    t_spr_perc::Array{Float64,1}
    t_spr_real::Array{Float64,1}
    t_eev::Array{Float64,1}
    t_ws::Array{Float64,1}
end