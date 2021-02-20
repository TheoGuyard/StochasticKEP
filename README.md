# Stochastic KEP

Code that is linked to the **[technical report](Stochastic_KEP.pdf)**.

## General informations

* Language : `Julia v1.4`
* Dependencies : `ArgParse`, `Cairo`, `Compose`,  `CSV`, `DataFrames`, `DelimitedFiles`, `Distributions`, `Fontconfig`, `GraphPlot`, `Gurobi`, `JLD2`, `JuMP`, `LightGraphs`, `MetaGraphs`, `Plots`, `PyPlot` 

## Usage

The file `exp/example.jl` contains examples on how to use `StochasticKep.jl`.

## Run a script

```bash
$ julia benchmarks/benchmark.jl "small.csv" [--nosave]
$ julia benchmarks/budget.jl "budget.csv" [--nosave]
```
