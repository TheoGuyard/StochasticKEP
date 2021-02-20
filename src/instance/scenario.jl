"""
    Scenario

General structure of a scenario for a KEP pool. The generation method for the 
oucomes is similar to the one in https://github.com/JohnDickerson/KidneyExchange.

# Parameters
* `pool::MetaDiGraph` : The KEP pool
* `generation::String` : Failure generation type
"""
mutable struct Scenario
    proba::Real
    generation::String
    failure_rate::Array{Float64,1}
    outcome::Array{Float64,1}
    function Scenario(proba::Real, generation::String, failure_rate::Array{Float64,1}, 
        outcome::Array{Float64,1})
        return new(proba, generation, failure_rate, outcome)
    end
    function Scenario(pool::MetaDiGraph, proba::Real, generation::String)
        0 <= proba <= 1 || throw(ArgumentError("Parameter proba must be in [0,1]."))
        generation in GENERATIONS || throw(ArgumentError("Generation not understood. Allowed generatations : $GENERATIONS."))
        failure_rate = get_failure_rate(pool, generation)
        outcome = zeros(ne(pool))
        for e in 1:ne(pool)
            f = failure_rate[e]
            outcome[e] = rand(Bernoulli(1-f)) ? 1. : 0.
        end
        return new(proba, generation, failure_rate, outcome)
    end
end

"""
    sample_scenarios

Generate a n-sample of scenarios with equal probability for given KEP pool and 
generation.

# Parameters
* `pool::MetaDiGraph` : The KEP pool
* `generation::String` : Failure generation type
* `n::Int` : Number of scenarios to generate
"""
function sample_scenarios(pool::MetaDiGraph, generation::String, n::Int)
    proba = 1/n
    scenarios = Scenario[]
    for s in 1:n
        push!(scenarios, Scenario(pool, proba, generation))
    end
    return scenarios
end

"""
    expected_scenario

Return the expected scenario of a sample of scenarios.

# Parameters
* `scenario::Array{Scenario,1}` : The sample of scenarios
"""
function expected_scenario(scenarios::Array{Scenario,1})
    length(scenarios) > 0 || throw(ArgumentError("Empty scenario array"))
    all(y -> y.generation == scenarios[1].generation, scenarios) || throw(ArgumentError("Scenarios do not have the same generation method"))
    mean_outcome = mean(scenario.outcome for scenario in scenarios)
    mean_failure = mean(scenario.failure_rate for scenario in scenarios)
    return Scenario(1., scenarios[1].generation, mean_failure, mean_outcome)
end

function get_failure_rate(pool::MetaDiGraph, generation::String)

    failure_rates = []

    for edge in edges(pool)
        # Dummies edges with an altruist as destination cannot fail
        if get_prop(pool, edge.dst, :altruist)
            push!(failure_rates, 0.)
        else
            if generation == "Constant"
                # Constant failure rate of 70%
                push!(failure_rates, 0.7)
            elseif generation == "Binomial"
                if rand() < 0.25
                    # Failure rate ~ 10%
                    push!(failure_rates, rand() * 0.2)
                else
                    # Failure rate ~ 90%
                    push!(failure_rates, 0.8 + rand() * 0.2)
                end
            elseif generation == "BinomialUNOS"
                # %pra of the patient < 0.8 : UNOS low sensitized patients
                if get_prop(pool, edge.dst, :pra) < 0.8
                    # Failure rate of 10%
                    push!(failure_rates, 0.1)
                else
                    # Failure rate of 90%
                    push!(failure_rates, 0.9)
                end
            elseif generation  == "BinomialAPD"
                # %pra of the patient < 0.75 : APD low sensitized patients
                if get_prop(pool, edge.dst, :pra) < 0.75
                    # Failure rate of 28%
                    push!(failure_rates, 0.28)
                else
                    # Failure rate of 58%
                    push!(failure_rates, 0.58)
                end
            elseif generation == "NoFailure"
                # No failures
                push!(failure_rates, 0.)
            end
        end
    end

    return failure_rates

end

Base.print(io::IO, scenario::Scenario) = println(io, "KEP scenario with probability $(scenario.proba), with $(scenario.generation) generation and with $(length(scenario.outcome)) outcomes")
Base.show(io::IO, scenario::Scenario) = print(io, "KEP scenario with probability $(scenario.proba)")
Base.print(io::IO, scenarios::Array{Scenario,1}) = println(io, "Array of $(length(scenarios)) KEP scenarios with $(scenarios[1].generation) generation and $(length(scenarios[1].outcome)) outcomes")
Base.show(io::IO, scenarios::Array{Scenario,1}) = print(io, "Array of $(length(scenarios)) KEP scenarios")