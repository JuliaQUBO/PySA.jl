# PySA.jl

[![CI](https://github.com/JuliaQUBO/PySA.jl/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/JuliaQUBO/PySA.jl/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/621844685.svg)](https://zenodo.org/badge/latestdoi/621844685)
[![QUBODRIVERS](https://img.shields.io/badge/Powered%20by-QUBODrivers.jl-%20%234063d8)](https://github.com/JuliaQUBO/QUBODrivers.jl)

[PySA](https://github.com/nasa/pysa) Simulated Annealing Interface for JuMP

## Installation
```julia
julia> import Pkg; Pkg.add("PySA")

julia> using PySA
```

## Getting started
```julia
using JuMP
using PySA

model = Model(PySA.Optimizer)
set_silent(model)
set_attribute(model, PySA.NumberOfReads(), 20)
set_attribute(model, PySA.NumberOfSweeps(), 64)

n = 3
Q = [ -1  2  2
       2 -1  2
       2  2 -1 ]

@variable(model, x[1:n], Bin)
@objective(model, Min, x' * Q * x)

optimize!(model)

for i = 1:result_count(model)
    xi = value.(x; result = i)
    yi = objective_value(model; result = i)
    println("[$i] f($(xi)) = $(yi)")
end
```

## Solver options

PySA.jl exposes the main PySA simulated annealing options as JuMP optimizer attributes:

| Attribute | PySA option | Default |
| --- | --- | --- |
| `PySA.NumberOfSweeps()` | `n_sweeps` | `32` |
| `PySA.NumberOfReplicas()` | `n_replicas` | `3` |
| `PySA.NumberOfReads()` | `n_reads` | `10` |
| `PySA.MinimumTemperature()` | `min_temp` | `1.0` |
| `PySA.MaximumTemperature()` | `max_temp` | `3.5` |
| `PySA.UpdateStrategy()` | `update_strategy` | `"sequential"` |
| `PySA.InitializeStrategy()` | `initialize_strategy` | `"ones"` |
| `PySA.RecomputeEnergy()` | `recompute_energy` | `false` |
| `PySA.SortOutputTemps()` | `sort_output_temps` | `true` |
| `PySA.Parallel()` | `parallel` | `true` |

Use `set_attribute(model, attribute, value)` to override these values before calling `optimize!`.
Use `set_silent(model)` to disable PySA solver output.

**Note**: _The PySA wrapper for Julia is not officially supported by the National Aeronautics and Space Administration. If you are interested in official support for Julia from NASA, let them know!_


**Note**: _If you are using `PySA.jl` in your project, we recommend you to include the `.CondaPkg` entry in your `.gitignore` file. The [`PythonCall`](https://github.com/cjdoris/PythonCall.jl) module will place a lot of files in this folder when building its Python environment._
