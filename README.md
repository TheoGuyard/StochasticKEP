# Stochastic KEP

Code that is linked to the **[technical report](Stochastic_KEP.pdf)**.

##### Abstract 
The Chronic Kidney Disease is the 11th most common cause of death globally, accounting for almost 1.2 million deaths worldwide per year. The most effective treatment is to receive a kidney transplant from another person. In many countries, kidney exchange schedules are done periodically. Pairs of donors and patients are pooled together, and the aim is to carry out as many transplants as possible. The problem is that in order to do a transplant, the donor and the recipient must be compatible. However, it is often not possible to test the compatibility between all the pairs in the pool. For this reason, strategies must be defined to choose which people to test in such a way that as many transplants as possible can be ultimately performed. These solutions result in long cycles of transplants involving several pairs, but a single failure due to incompatibility on one transplant causes the entire cycle in which it is included to fail. It is thus necessary to take into account the fact that all the transplants will not necessarily be able to be done but we do not know a priori which ones will be. In this project, stochastic programming models are used to tackle the problem. These models are designed to take into account the uncertainty of the data in the problem. In particular, two models corresponding to different practical cases are derived. The quality of these two models is evaluated along with the gains compared to other strategies, for instance the one where we only consider the average compatibility failure to decide about the strategy. The numerical results that are obtained will allow us to define strategies according to the number of tests that can be performed in the pool.

## General informations

* Language : `Julia v1.4`
* Dependencies : `ArgParse`, `DataFrames`, `Distributions`, `GLPK`, `JuMP`, `LightGraphs`, `MetaGraphs`, `Plots`

## Install packages

To install the required packages (in a virtual environnement), open a command line terminal in the folder `StochasticKEP` and run the following commands in Julia's REPL.
```julia
julia> ]
(@v1.4) pkg> activate .
(@v1.4) pkg> instantiate
```
Models are solved using `GLPK`.

## Data files

Data files are located in the `StochasticKEP/data` folder. They come from the [PrefLib](https://www.preflib.org) library. In this repository, only one dataset in included. Additional datasets can be downloaded from [here](https://www.preflib.org/data/matching/kidney/) and added manually in the data folder.

## Example scripts

The following commands can be run from the `StochasticKEP` folder.

### Allowed command-line arguments values

* `instance` : relative path of the data set from `StochasticKEP/data` folder
* `formulation` : `matching`, `relaxed_arc`
* `generation` : `Constant`, `Binomial`, `BinomialUNOS`, `BinomialAPD`, `NoFailure`
* `scenarios` : Positive `Int`
* `eval-scenarios` : Positive `Int`

### Collect information about a dataset
Print the graph size for the matching and the relaxed-arc formulation with the corresponding number of variables in the models.
```bash
$ julia --project=. example/instance.jl <instance> <scenarios>
$ # Example : julia --project=. example/instance.jl preflib-md/MD-00001-00000001 100
```

### Perform a scenario sensibility analysis
For different numbers of scenarios in a given formulation, plot the real cost, the perceived cost and the solution time.
```bash
$ julia --project=. example/scenarios.jl <instance> <formulation> <generation> <budget> <min-scenario> <max-scenario> <step-scenario> <eval-scenarios> <repeats> [<maxtime>]
$ # Example : julia --project=. example/scenarios.jl preflib-md/MD-00001-00000001 matching BinomialUNOS 4 1 10 2 100 3 60
```
Results are stored in `example/saves/scenario.png`.

### Perform a budget sensibility analysis
For different HLA-crossmatch test budgets in a given formulation, plot the real cost, the perceived cost, the EEV, the WS and the solution time.
```bash
$ julia --project=. example/budget.jl <instance> <formulation> <generation> <min-budget> <max-budget> <step-budget> <scenarios> <eval-scenarios> [<maxtime>]
$ # Example : julia --project=. example/budget.jl preflib-md/MD-00001-00000001 matching BinomialUNOS 1 10 2 10 100 60
```
Results are stored in `example/saves/budget.png`.
