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

println("\n------")
println("Budget")
println("------")

for i in 1:nb_instances

    try
        # Extract instance specifications
        instance = setupfile[i,:instance]
        formulation = setupfile[i,:formulation]
        generation = setupfile[i,:generation]
        K = setupfile[i,:K]
        K_eval = setupfile[i,:K_eval]
        min_budget_factor = setupfile[i,:min_bf]
        max_budget_factor = setupfile[i,:max_bf]
        nb_budget_factors = setupfile[i,:nb_bf]
        maxtime = setupfile[i,:maxtime]

        step = (max_budget_factor-min_budget_factor)/(nb_budget_factors-1)
        budget_factors = min_budget_factor:step:max_budget_factor

        # Results
        z_sp_perc = []
        z_sp_real = []
        z_spr_perc = []
        z_spr_real = []
        eev = []
        ws = []
        t_sp_perc = []
        t_sp_real = []
        t_spr_perc = []
        t_spr_real = []
        t_eev = []
        t_ws = []

        # Construct instance
        wmd_file = joinpath(data_folder, join([instance, ".wmd"]))
        dat_file = joinpath(data_folder, join([instance, ".dat"]))
        pool = kep_pool(wmd_file, dat_file)

        # Sample scenarios
        scenarios = sample_scenarios(pool, generation, K)
        eval_scenarios = sample_scenarios(pool, generation, K_eval)

        println("\nTreating instance $i/$(nb_instances) ...")
        println("  instance\t: $instance")
        println("  formulation\t: $formulation")
        println("  generation\t: $generation")
        println("  scenarios\t: $(K)")
        println("  ev. scenarios\t: $(K_eval)")
        println("  maxtime\t: $maxtime seconds")
        println()

        # Reference solution for budget factors
        ref_scenario = Scenario(pool, 1, "NoFailure")
        ref_var, ref_monitor = DP(formulation, pool, ref_scenario, ne(pool), maxtime=maxtime)
        z_max = ref_monitor["objective_value"]
        z_max != -Inf || println("    \e[91mComputation of z_max exceed maxtime\e[00m")

        for (b, budget_factor) in enumerate(budget_factors)

            println("  Budget factor : $(budget_factor) ...")
        
            budget = convert(Int, ceil(budget_factor * z_max))
        
            v_perc, m_perc = SP(formulation, pool, scenarios, budget, maxtime=maxtime)
            v_real, m_real = decision_evaluation(formulation, v_perc["x"], pool, eval_scenarios, maxtime=maxtime)
            vr_perc, mr_perc = SP(formulation, pool, scenarios, budget, maxtime=maxtime, relax=true)
            vr_real, mr_real = decision_evaluation(formulation, vr_perc["x"], pool, eval_scenarios, maxtime=maxtime, relax=true)
            v_evp, m_evp = EVP(formulation, pool, scenarios, budget, maxtime=maxtime)
            v_eev, m_eev = decision_evaluation(formulation, v_evp["x"], pool, eval_scenarios, maxtime=maxtime)
            wsv, m_ws = WS(formulation, pool, eval_scenarios, budget, maxtime=maxtime)

            push!(z_sp_perc, m_perc["objective_value"])
            push!(z_sp_real, m_real["objective_value"])
            push!(z_spr_perc, mr_perc["objective_value"])
            push!(z_spr_real, mr_real["objective_value"])
            push!(eev, m_eev["objective_value"])
            push!(ws, wsv)
            push!(t_sp_perc, m_perc["solve_time"])
            push!(t_sp_real, m_real["solve_time"])
            push!(t_spr_perc, mr_perc["solve_time"])
            push!(t_spr_real, mr_real["solve_time"])
            push!(t_eev, m_eev["solve_time"])
            push!(t_ws, m_ws["solve_time"])
        end

        budget_instance = BudgetInstance(
            pool,
            formulation,
            generation,
            K,
            K_eval,
            scenarios,
            eval_scenarios,
            maxtime,
            budget_factors,
            z_sp_perc,
            z_sp_real,
            z_spr_perc,
            z_spr_real,
            eev,
            ws,
            t_sp_perc,
            t_sp_real,
            t_spr_perc,
            t_spr_real,
            t_eev,
            t_ws
        )

        title = join([
            "benchmarks/results/",
            "budget_",
            split(instance, "/")[end],
            "_",
            generation,
            "_",
            formulation,
            "_",
            Dates.format(now(), "yyyy-mm-dd_HH-MM-SS"),
            ".jld2"
        ])

        if !parsed_args["nosave"]
            @save title budget_instance
        end
    catch e
        println("\e[91mError with instance $i\e[00m")
        println(e)
        continue
    end
end