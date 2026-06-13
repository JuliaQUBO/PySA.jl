module PySA

using PythonCall
using LinearAlgebra

import QUBODrivers
import QUBOTools
import MathOptInterface as MOI

const np      = PythonCall.pynew()
const pysa    = PythonCall.pynew()
const pysa_sa = PythonCall.pynew()
const seed_numba_rng = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(pysa, pyimport("pysa"))
    PythonCall.pycopy!(pysa_sa, pyimport("pysa.sa"))

    ns = PythonCall.pydict()
    pyimport("builtins").exec(
        """
        import numpy as np
        import numba

        @numba.njit
        def seed_numba_rng(seed):
            np.random.seed(seed)
        """,
        ns,
    )
    PythonCall.pycopy!(seed_numba_rng, ns["seed_numba_rng"])
end

QUBODrivers.@setup Optimizer begin
    name       = "PySA"
    version    = v"0.1.0" # pysa version
    attributes = begin
        NumberOfSweeps["num_sweeps"]::Integer = 32
        NumberOfReplicas["num_replicas"]::Integer = 3
        NumberOfReads["num_reads"]::Integer = 10
        "seed"::Union{Integer,Nothing} = nothing
        MinimumTemperature["min_temp"]::Float64 = 1.0
        MaximumTemperature["max_temp"]::Float64 = 3.5
        UpdateStrategy["update_strategy"]::String = "sequential"
        InitializeStrategy["initialize_strategy"]::String = "ones"
        RecomputeEnergy["recompute_energy"]::Bool = false
        SortOutputTemps["sort_output_temps"]::Bool = true
        Parallel["parallel"]::Bool = true
    end
end

const RandomSeed = QUBODrivers.RandomSeed
const FinalNumberOfReads = QUBODrivers.FinalNumberOfReads

const _RawLegacyNumberOfSweeps = QUBODrivers.RawSamplerAttribute{Symbol("n_sweeps")}
const _RawLegacyNumberOfReplicas = QUBODrivers.RawSamplerAttribute{Symbol("n_replicas")}
const _RawLegacyNumberOfReads = QUBODrivers.RawSamplerAttribute{Symbol("n_reads")}
const _RawMinimumTemperature = QUBODrivers.RawSamplerAttribute{Symbol("minimum_temperature")}
const _RawMaximumTemperature = QUBODrivers.RawSamplerAttribute{Symbol("maximum_temperature")}

MOI.get(sampler::Optimizer, ::_RawLegacyNumberOfSweeps) = MOI.get(sampler, NumberOfSweeps())
MOI.get(sampler::Optimizer, ::_RawLegacyNumberOfReplicas) = MOI.get(sampler, NumberOfReplicas())
MOI.get(sampler::Optimizer, ::_RawLegacyNumberOfReads) = MOI.get(sampler, NumberOfReads())
MOI.get(sampler::Optimizer, ::_RawMinimumTemperature) = MOI.get(sampler, MinimumTemperature())
MOI.get(sampler::Optimizer, ::_RawMaximumTemperature) = MOI.get(sampler, MaximumTemperature())

function MOI.set(sampler::Optimizer, ::_RawLegacyNumberOfSweeps, value)
    MOI.set(sampler, NumberOfSweeps(), value)
    return nothing
end

function MOI.set(sampler::Optimizer, ::_RawLegacyNumberOfReplicas, value)
    MOI.set(sampler, NumberOfReplicas(), value)
    return nothing
end

function MOI.set(sampler::Optimizer, ::_RawLegacyNumberOfReads, value)
    MOI.set(sampler, NumberOfReads(), value)
    return nothing
end

function MOI.set(sampler::Optimizer, ::_RawMinimumTemperature, value)
    MOI.set(sampler, MinimumTemperature(), value)
    return nothing
end

function MOI.set(sampler::Optimizer, ::_RawMaximumTemperature, value)
    MOI.set(sampler, MaximumTemperature(), value)
    return nothing
end

MOI.supports(::Optimizer, ::_RawLegacyNumberOfSweeps) = true
MOI.supports(::Optimizer, ::_RawLegacyNumberOfReplicas) = true
MOI.supports(::Optimizer, ::_RawLegacyNumberOfReads) = true
MOI.supports(::Optimizer, ::_RawMinimumTemperature) = true
MOI.supports(::Optimizer, ::_RawMaximumTemperature) = true

QUBODrivers.honors_final_reads(::Type{<:Optimizer}) = true
QUBODrivers.enforces_time_limit(::Type{<:Optimizer}) = false

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

function _check_nonnegative_integer(name::String, value::Integer)
    value >= 0 || error("Value for '$name' must be a non-negative integer")

    return nothing
end

function _seed_backend!(::Nothing)
    return nothing
end

function _seed_backend!(seed::Integer)
    np.random.seed(seed)
    seed_numba_rng(seed)

    return nothing
end

function _pysa_spin_sample(state, h, J, α, β, ::Type{T}) where {T}
    ψ = -round.(Int, pyconvert.(T, [var for var in state]))
    λ = QUBOTools.value(ψ, h, J, α, β)

    return QUBOTools.Sample{T,Int}(ψ, λ)
end

function _final_samples(result, final_num_reads::Integer, h, J, α, β, ::Type{T}) where {T}
    best_states = result["best_state"].values
    backend_num_reads = length(best_states)
    samples = Vector{QUBOTools.Sample{T,Int}}(undef, backend_num_reads)

    for i = 1:backend_num_reads
        # NOTE: Python is 0-indexed
        # NOTE: sign has to be inverted, as mentioned before
        samples[i] = _pysa_spin_sample(best_states[i-1], h, J, α, β, T)
    end

    sort!(samples; by = QUBOTools.value)

    return samples[1:final_num_reads]
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

    num_sweeps = MOI.get(sampler, PySA.NumberOfSweeps())
    num_replicas = MOI.get(sampler, PySA.NumberOfReplicas())
    num_reads = MOI.get(sampler, PySA.NumberOfReads())
    final_num_reads = MOI.get(sampler, QUBODrivers.FinalNumberOfReads())
    seed = MOI.get(sampler, QUBODrivers.RandomSeed())
    requested_parallel = MOI.get(sampler, PySA.Parallel())
    backend_parallel = requested_parallel && isnothing(seed)
    backend_num_reads = max(num_reads, final_num_reads)

    _check_nonnegative_integer("NumberOfSweeps", num_sweeps)
    _check_nonnegative_integer("NumberOfReplicas", num_replicas)
    _check_nonnegative_integer("NumberOfReads", num_reads)
    _check_nonnegative_integer("FinalNumberOfReads", final_num_reads)
    num_replicas > 0 || error("Value for 'NumberOfReplicas' must be positive")

    result = @timed begin
        if backend_num_reads == 0
            nothing
        else
            _seed_backend!(seed)
            solver.metropolis_update(
                num_sweeps          = num_sweeps,
                num_reads           = backend_num_reads,
                num_replicas        = num_replicas,
                update_strategy     = MOI.get(sampler, PySA.UpdateStrategy()),
                min_temp            = MOI.get(sampler, PySA.MinimumTemperature()),
                max_temp            = MOI.get(sampler, PySA.MaximumTemperature()),
                initialize_strategy = MOI.get(sampler, PySA.InitializeStrategy()),
                recompute_energy    = MOI.get(sampler, PySA.RecomputeEnergy()),
                sort_output_temps   = MOI.get(sampler, PySA.SortOutputTemps()),
                parallel            = backend_parallel,
                verbose             = !MOI.get(sampler, MOI.Silent()),
            )
        end
    end

    samples = if backend_num_reads == 0
        QUBOTools.Sample{T,Int}[]
    else
        _final_samples(result.value, final_num_reads, h, J, α, β, T)
    end

    seeds = if isnothing(seed)
        Dict{String,Any}()
    else
        Dict{String,Any}("numpy" => seed, "numba" => seed)
    end

    metadata = QUBODrivers._sampler_metadata(
        origin                = "PySA.jl",
        algorithm_name        = "simulated_annealing",
        backend_name          = "pysa",
        backend_version       = MOI.get(sampler, MOI.SolverVersion()),
        execution_mode        = backend_parallel ? "parallel" : "sequential",
        optimizer_iterations  = num_sweeps,
        optimizer_evaluations = backend_num_reads * num_replicas,
        number_of_reads       = backend_num_reads,
        final_number_of_reads = final_num_reads,
        seeds                 = seeds,
        status                = "locally_solved",
        termination_status    = MOI.LOCALLY_SOLVED,
    )
    metadata["time"] = Dict{String,Any}("effective" => result.time)
    metadata["parameters"] = Dict{String,Any}(
        "num_sweeps"         => num_sweeps,
        "num_replicas"       => num_replicas,
        "num_reads"          => num_reads,
        "final_num_reads"    => final_num_reads,
        "requested_parallel" => requested_parallel,
        "parallel"           => backend_parallel,
    )

    return QUBOTools.SampleSet{T}(samples, metadata; domain = :spin, sense = :min)
end

end # module PySA
