using ArgParse, CSV, Dates, DataFrames, JLD2, JuMP, LightGraphs, MetaGraphs

include("../src/StochasticKep.jl")
using .StochasticKep

include("typedef.jl")

data_folder = joinpath(join(split(@__DIR__, "/")[1:(end-1)], "/"), "data")


# ----- Argument parsing ----- #

argparse_settings = ArgParseSettings()
@add_arg_table! argparse_settings begin
    "setupfile"
        help = "the setup file of the benchmarks"
        arg_type = String
        required = true
    "--nosave"
        help = "not saving results"
        action = :store_true
end

parsed_args = parse_args(ARGS, argparse_settings)


# ----- Benchmark ----- #

setupfile = DataFrame(CSV.File(join(["benchmarks/setupfiles/", parsed_args["setupfile"]])))
nb_instances = size(setupfile)[1]

println("\n----------")
println("Benchmarks")
println("----------")

for i in 1:nb_instances

    try
        # Extract instance specifications
        instance = setupfile[i,:instance]
        formulation = setupfile[i,:formulation]
        budget = setupfile[i,:budget]
        generation = setupfile[i,:generation]
        K = setupfile[i,:K]
        K_eval = setupfile[i,:K_eval]
        maxtime = setupfile[i,:maxtime]

        maxcyclelength = 4

        # Construct instance
        wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
        dat_file = joinpath(data_folder, join([instance, ".dat"]))
        pool = kep_pool(wmd_file, dat_file)

        # Generate scenarios
        scenarios = sample_scenarios(pool, generation, K)
        eval_scenarios = sample_scenarios(pool, generation, K_eval)

        println("\nTreating instance $i/$(nb_instances) ...")
        println("  instance\t: $instance")
        println("  formulation\t: $formulation")
        println("  budget\t: $budget")
        println("  generation\t: $generation")
        println("  scenarios\t: $(K)")
        println("  ev. scenarios\t: $(K_eval)")
        println("  maxtime\t: $maxtime seconds")
        println()

        results = Dict{String,Any}()
        
        println("  Solving the stochastic problem ...")
        v_sp, m_sp = SP(formulation, pool, scenarios, budget, maxtime=maxtime)
        v_sp_eval, m_sp_eval = decision_evaluation(formulation, v_sp["x"], pool, eval_scenarios, maxtime=maxtime)
        results["sp"] = Dict(
            "vars" => v_sp,
            "monitoring" => m_sp,
            "perc_cost" => m_sp["objective_value"],
            "real_cost" => m_sp_eval["objective_value"],
        )

        println("  Finding cycles ...")
        cycles, mean_nb_cycle = find_cycles(v_sp["y"], pool, maxlength=maxcyclelength)
        results["mean_nb_cycles"] = mean_nb_cycle

        println("  Solving the stochastic problem (relaxed) ...")
        v_spr, m_spr = SP(formulation, pool, scenarios, budget, maxtime=maxtime, relax=true)
        v_spr_eval, m_spr_eval = decision_evaluation(formulation, v_spr["x"], pool, eval_scenarios, maxtime=maxtime)
        results["spr"] = Dict(
            "vars" => v_spr,
            "monitoring" => m_spr,
            "perc_cost" => m_spr["objective_value"],
            "real_cost" => m_spr_eval["objective_value"],
        )

        println("  Computing the EEV ...")
        v_evp, m_evp = EVP(formulation, pool, scenarios, budget, maxtime=maxtime)
        v_eev, m_eev = decision_evaluation(formulation, v_evp["x"], pool, eval_scenarios, maxtime=maxtime)
        results["eev"] = Dict(
            "vars" => v_evp,
            "monitoring" => m_evp,
            "perc_cost" => m_evp["objective_value"],
            "real_cost" => m_eev["objective_value"],
        )

        println("  Computing the WS ...")
        ws, m_ws = WS(formulation, pool, eval_scenarios, budget, maxtime=maxtime)
        results["ws"] = Dict(
            "real_cost" => ws,
            "monitoring" => m_ws,
        )

        benchmark_instance = BenchmarkInstance(
            pool,
            formulation,
            budget,
            generation, 
            K,
            K_eval,
            scenarios,
            eval_scenarios,
            maxtime,
            maxcyclelength,
            results
        )

        title = join([
            "benchmarks/results/",
            split(instance, "/")[end],
            "_",
            formulation,
            "_",
            Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"),
            ".jld2"
        ])

        if !parsed_args["nosave"]
            @save title benchmark_instance
        end
    catch e
        println("\e[91mError with instance $i\e[00m")
        println(e)
        continue
    end
end





