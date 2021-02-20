using CSV, DataFrames, Plots

file = "2021-01-23_13-00-24.csv"

table = DataFrame(CSV.File(join(["benchmarks/saves/", file])))

table_m = table[table.formulation .== "matching", :]
table_r = table[table.formulation .== "relaxed_arc", :]
table_m_1 = table_m[table_m.budget .== table_m.V, :]
table_m_2 = table_m[table_m.budget .== 2*table_m.V, :]
table_r_1 = table_r[table_r.budget .== table_r.V, :]
table_r_2 = table_r[table_r.budget .== 2*table_r.V, :]

instances_m_1 = unique(table_m_1.instance)
instances_m_2 = unique(table_m_2.instance)
instances_r_1 = unique(table_r_1.instance)
instances_r_2 = unique(table_r_2.instance)

size_m_1 = 1:length(instances_m_1)
size_m_2 = 1:length(instances_m_2)
size_r_1 = 1:length(instances_r_1)
size_r_2 = 1:length(instances_r_2)

# --- Quality --- #

plot(
    size_m_2,
    [table_m_2.z_sp_perc, table_m_2.z_sp_real, table_m_2.z_eev_real, table_m_2.ws, table_m_2.z_spr_real],
    label = ["perc cost" "real cost" "EEVr" "WS" "real cost (relax)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 6,
    legend = :topleft,
    legendfontsize = 11,
    xlabel = "Instance index",
    ylabel = "Objective",
    xticks = size_m_2,
    guidefontsize = 13,
    size = (600, 400)
)

# --- Solve times --- #

plot(
    size_r_2,
    [table_r_2.t_sp, table_r_2.t_spr],
    label = ["t" "t (relax)"],
    yaxis = :log,
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 6,
    legend = :topleft,
    legendfontsize = 11,
    xlabel = "Instance index",
    xticks = size_r_2,
    ylabel = "Seconds",
    guidefontsize = 13,
    size = (600, 400)
)

# --- EVPI --- #

plot(
    [size_m_1, size_r_1],
    [table_m_1.ws - table_m_1.z_sp_real, table_r_1.ws - table_r_1.z_sp_real],
    label = ["EVPI - B=|V| (matching)" "EVPI - B=|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :topright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (1200, 300),
    left_margin = 5Plots.mm, 
    bottom_margin = 5Plots.mm
)

plot!(
    [size_m_2, size_r_2],
    [table_m_2.ws - table_m_2.z_sp_real, table_r_2.ws - table_r_2.z_sp_real],
    label = ["EVPI - B=2|V| (matching)" "EVPI - B=2|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :utriangle,
    markersize = 5,
    legend = :topright,
)


# --- VSS --- #

plot(
    [size_m_1, size_r_1],
    [table_m_1.z_sp_real - table_m_1.z_eev_real, table_r_1.z_sp_real - table_r_1.z_eev_real],
    label = ["VSS - B=|V| (matching)" "VSS - B=|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :topright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (1200, 300),
    left_margin = 5Plots.mm, 
    bottom_margin = 5Plots.mm
)

plot!(
    [size_m_2, size_r_2],
    [table_m_2.z_sp_real - table_m_2.z_eev_real, table_r_2.z_sp_real - table_r_2.z_eev_real],
    label = ["VSS - B=2|V| (matching)" "VSS - B=2|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :utriangle,
    markersize = 5,
    legend = :topright,
)

# --- VRS (Value of the Relaxed Solution) --- #

plot(
    [size_m_1, size_r_1],
    [table_m_1.z_spr_real - table_m_1.z_sp_real, table_r_1.z_spr_real - table_r_1.z_sp_real],
    label = ["VRS - B=|V| (matching)" "VRS - B=|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :topright,
    xlabel = "Instance index",
    ylabel = "Value",
    size = (1200, 300),
    left_margin = 5Plots.mm, 
    bottom_margin = 5Plots.mm
)

plot!(
    [size_m_2, size_r_2],
    [table_m_2.z_spr_real - table_m_2.z_sp_real, table_r_2.z_spr_real - table_r_2.z_sp_real],
    label = ["VSS - B=2|V| (matching)" "VSS - B=2|V| (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :utriangle,
    markersize = 5,
    legend = :topright,
)


# --- Cycles --- #

p1 = plot(
    [size_m_1, size_r_1, size_r_1, size_r_1],
    [table_m_1.cycle_2, table_r_1.cycle_2, table_r_1.cycle_3, table_r_1.cycle_4_plus],
    label = ["C_2 (matching)" "C_2 (relaxed-arc)" "C_3 (relaxed-arc)" "C_4+ (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    title = "Budget = |V|",
    xlabel = "Instance index",
    ylabel = "Mean number of cycles",
    size = (1200, 300),
    left_margin = 5Plots.mm, 
    bottom_margin = 5Plots.mm
)

p2 = plot(
    [size_m_2, size_r_2, size_r_2, size_r_2],
    [table_m_2.cycle_2, table_r_2.cycle_2, table_r_2.cycle_3, table_r_2.cycle_4_plus],
    label = ["C_2 (matching)" "C_2 (relaxed-arc)" "C_3 (relaxed-arc)" "C_4+ (relaxed-arc)"],
    lw = 2,
    ls = :dash,
    markershape = :circle,
    markersize = 5,
    legend = :bottomright,
    title = "Budget = 2|V|",
    xlabel = "Instance index",
    ylabel = "Mean number of cycles",
    size = (1200, 300),
    left_margin = 5Plots.mm, 
    bottom_margin = 5Plots.mm
)

plot(p1, p2, layout=(2,1), size=(1200,600))
