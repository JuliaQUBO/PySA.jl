# PySA.jl

[![DOI](https://zenodo.org/badge/621844685.svg)](https://zenodo.org/badge/latestdoi/621844685)
[![QUBODRIVERS](https://img.shields.io/badge/Powered%20by-QUBODrivers.jl-%20%234063d8)](https://github.com/psrenergy/QUBODrivers.jl)

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

**Note**: _The PySA wrapper for Julia is not officially supported by the National Aeronautics and Space Administration. If you are interested in official support for Julia from NASA, let them know!_


**Note**: _If you are using `PySA.jl` in your project, we recommend you to include the `.CondaPkg` entry in your `.gitignore` file. The [`PythonCall`](https://github.com/cjdoris/PythonCall.jl) module will place a lot of files in this folder when building its Python environment._
