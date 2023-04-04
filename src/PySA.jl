module PySA

using PythonCall
using LinearAlgebra
using Anneal

const np = PythonCall.pynew()
const pysa = PythonCall.pynew()
const pysa_sa = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(np, pyimport("numpy"))
    PythonCall.pycopy!(pysa, pyimport("pysa"))
    PythonCall.pycopy!(pysa_sa, pyimport("pysa.sa"))
end

Anneal.@anew Optimizer begin
    name       = "PySA"
    sense      = :min
    domain     = :spin
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

function Anneal.sample(sampler::Optimizer{T}) where {T}
    h, J, α, β = Anneal.ising(sampler, Matrix)

    problem = np.array(Symmetric(J + diagm(h)))

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

    samples = Anneal.Sample{T,Int}[]

    for (Ψ, Λ) in result.value[pylist(["states", "energies"])].values
        for i = 0:(num_replicas-1)
            ψ = round.(Int, pyconvert.(T, Ψ[i]))
            λ = α * (pyconvert(T, Λ[i]) + β)

            push!(samples, Anneal.Sample{T}(ψ, λ))
        end
    end

    metadata = Dict{String,Any}(
        "time" => Dict{String,Any}(
            "effective" => result.time
        )
    )

    return Anneal.SampleSet{T}(samples, metadata)
end

end # module PySA
