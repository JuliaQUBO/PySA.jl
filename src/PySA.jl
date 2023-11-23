module PySA

using PythonCall
using LinearAlgebra

import QUBODrivers
import QUBOTools
import MathOptInterface as MOI

const np      = PythonCall.pynew()
const pysa    = PythonCall.pynew()
const pysa_sa = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(pysa, pyimport("pysa"))
    PythonCall.pycopy!(pysa_sa, pyimport("pysa.sa"))
end

QUBODrivers.@setup Optimizer begin
    name       = "PySA"
    version    = v"0.1.0" # pysa version
    attributes = begin
        NumberOfSweeps["n_sweeps"]::Integer = 32
        NumberOfReplicas["n_replicas"]::Integer = 3
        NumberOfReads["n_reads"]::Integer = 10
        MinimumTemperature["min_temp"]::Float64 = 1.0
        MaximumTemperature["max_temp"]::Float64 = 3.5
        UpdateStrategy["update_strategy"]::String = "sequential"
        InitializeStrategy["initialize_strategy"]::String = "ones"
        RecomputeEnergy["recompute_energy"]::Bool = false
        SortOutputTemps["sort_output_temps"]::Bool = true
        Parallel["parallel"]::Bool = true
    end
end

function _float_type(::Type{T})::String where {T<:AbstractFloat}
    if T === Float16
        return "float16"
    elseif T === Float32
        return "float32"
    elseif T === Float64
        return "float64"
    else
        error("Unknown float type '$T'")
    end
end

function QUBODrivers.sample(sampler::Optimizer{T}) where {T}
    n, h, J, α, β = QUBOTools.ising(sampler, :dense; sense = :min)

    # Since PySA adopts the s = 1 - 2x instead of the s = 2x - 1
    # convention, the sign of 'h' has to be inverted, as well as
    # the value for each state variable in 'ψ' below
    problem = np.array(Symmetric(J - diagm(h)))

    solver = pysa_sa.Solver(
        problem      = problem,
        problem_type = "ising",
        float_type   = _float_type(T),
    )

    num_replicas = MOI.get(sampler, PySA.NumberOfReplicas())

    result = @timed solver.metropolis_update(
        num_sweeps          = MOI.get(sampler, PySA.NumberOfSweeps()),
        num_reads           = MOI.get(sampler, PySA.NumberOfReads()),
        num_replicas        = num_replicas,
        update_strategy     = MOI.get(sampler, PySA.UpdateStrategy()),
        min_temp            = MOI.get(sampler, PySA.MinimumTemperature()),
        max_temp            = MOI.get(sampler, PySA.MaximumTemperature()),
        initialize_strategy = MOI.get(sampler, PySA.InitializeStrategy()),
        recompute_energy    = MOI.get(sampler, PySA.RecomputeEnergy()),
        sort_output_temps   = MOI.get(sampler, PySA.SortOutputTemps()),
        parallel            = MOI.get(sampler, PySA.Parallel()), # True by default
        verbose             = !MOI.get(sampler, MOI.Silent()),
    )

    samples = Vector{QUBOTools.Sample{T,Int}}(undef, num_replicas)

    for Ψ in result.value["states"].values
        for i = 1:num_replicas
            # NOTE: Python is 0-indexed
            # NOTE: sign has to be inverted, as mentioned before
            ψ = -round.(Int, pyconvert.(T, Ψ[i-1]))

            # Compute value instead of retrieving it, to avoid precision errors
            λ = QUBOTools.value(ψ, h, J, α, β) 

            samples[i] = QUBOTools.Sample{T,Int}(ψ, λ)
        end
    end

    metadata = Dict{String,Any}(
        "origin" => "PySA",
        "time"   => Dict{String,Any}(
            "effective" => result.time
        ),
    )

    return QUBOTools.SampleSet{T}(samples, metadata; domain = :spin, sense = :min)
end

end # module PySA
