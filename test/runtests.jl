import PySA
import PySA: MOI, QUBODrivers, QUBOTools
import JuMP
import TOML
using Test

@testset "Package metadata" begin
    project = TOML.parsefile(joinpath(pkgdir(PySA), "Project.toml"))
    condapkg = TOML.parsefile(joinpath(pkgdir(PySA), "CondaPkg.toml"))
    qubodrivers_compat = split(project["compat"]["QUBODrivers"], r",\s*")
    readme = read(joinpath(pkgdir(PySA), "README.md"), String)
    source = read(joinpath(pkgdir(PySA), "src", "PySA.jl"), String)
    pysa_pin = match(r"^@ git\+https://github\.com/nasa/pysa@v(\d+\.\d+\.\d+)$", condapkg["pip"]["deps"]["pysa"])
    pysa_version = match(r"(?m)^\s*version\s*=\s*v\"(\d+\.\d+\.\d+)\"\s*(?:#.*)?$", source)

    @test project["compat"]["QUBODrivers"] == "0.6.1"
    @test project["compat"]["QUBOTools"] == "0.13, 0.14, 0.15, 0.16"
    @test !("0.5" in qubodrivers_compat)
    @test occursin("https://github.com/JuliaQUBO/QUBODrivers.jl", readme)
    @test !occursin("https://github.com/psrenergy/QUBODrivers.jl", readme)
    @test occursin("## Solver options", readme)
    for attribute in (
        "PySA.NumberOfSweeps()",
        "PySA.NumberOfReplicas()",
        "PySA.NumberOfReads()",
        "PySA.FinalNumberOfReads()",
        "PySA.RandomSeed()",
        "PySA.MinimumTemperature()",
        "PySA.MaximumTemperature()",
        "PySA.UpdateStrategy()",
        "PySA.InitializeStrategy()",
        "PySA.RecomputeEnergy()",
        "PySA.SortOutputTemps()",
        "PySA.Parallel()",
    )
        @test occursin(attribute, readme)
    end
    @test pysa_pin !== nothing
    @test pysa_version !== nothing

    if pysa_pin !== nothing
        @test pysa_pin.captures[1] == "0.1.0"
    end

    if pysa_pin !== nothing && pysa_version !== nothing
        @test pysa_version.captures[1] == pysa_pin.captures[1]
    end
end

@testset "Solver attributes" begin
    sampler = PySA.Optimizer()

    @test MOI.supports(sampler, PySA.RandomSeed())
    @test MOI.supports(sampler, PySA.FinalNumberOfReads())
    @test QUBODrivers.supports_seed(PySA.Optimizer)
    @test QUBODrivers.honors_final_reads(PySA.Optimizer)
    @test !QUBODrivers.enforces_time_limit(PySA.Optimizer)

    MOI.set(sampler, PySA.RandomSeed(), 123)
    @test MOI.get(sampler, PySA.RandomSeed()) == 123
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("seed")) == 123

    MOI.set(sampler, PySA.FinalNumberOfReads(), 4)
    @test MOI.get(sampler, PySA.FinalNumberOfReads()) == 4
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("final_num_reads")) == 4

    MOI.set(sampler, MOI.RawOptimizerAttribute("n_reads"), 3)
    @test MOI.get(sampler, PySA.NumberOfReads()) == 3
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("num_reads")) == 3

    MOI.set(sampler, MOI.RawOptimizerAttribute("num_reads"), 5)
    @test MOI.get(sampler, PySA.NumberOfReads()) == 5
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("n_reads")) == 5

    MOI.set(sampler, MOI.RawOptimizerAttribute("n_sweeps"), 7)
    @test MOI.get(sampler, PySA.NumberOfSweeps()) == 7
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("num_sweeps")) == 7

    MOI.set(sampler, MOI.RawOptimizerAttribute("n_replicas"), 2)
    @test MOI.get(sampler, PySA.NumberOfReplicas()) == 2
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("num_replicas")) == 2

    MOI.set(sampler, MOI.RawOptimizerAttribute("minimum_temperature"), 0.25)
    @test MOI.get(sampler, PySA.MinimumTemperature()) == 0.25
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("min_temp")) == 0.25

    MOI.set(sampler, MOI.RawOptimizerAttribute("maximum_temperature"), 4.0)
    @test MOI.get(sampler, PySA.MaximumTemperature()) == 4.0
    @test MOI.get(sampler, MOI.RawOptimizerAttribute("max_temp")) == 4.0
end

@testset "README JuMP workflow" begin
    model = JuMP.Model(PySA.Optimizer)
    JuMP.set_attribute(model, MOI.Silent(), true)
    JuMP.set_attribute(model, PySA.RandomSeed(), 123)
    JuMP.set_attribute(model, PySA.NumberOfReads(), 4)
    JuMP.set_attribute(model, PySA.FinalNumberOfReads(), 4)
    JuMP.set_attribute(model, PySA.NumberOfSweeps(), 8)
    JuMP.set_attribute(model, PySA.NumberOfReplicas(), 2)
    JuMP.set_attribute(model, PySA.MinimumTemperature(), 0.5)
    JuMP.set_attribute(model, PySA.MaximumTemperature(), 2.5)

    n = 3
    Q = [
        -1  2  2
         2 -1  2
         2  2 -1
    ]

    JuMP.@variable(model, x[1:n], Bin)
    JuMP.@objective(model, Min, x' * Q * x)

    JuMP.optimize!(model)

    @test JuMP.result_count(model) >= 1
    objective = JuMP.objective_value(model; result = 1)
    @test isfinite(objective)
    @test objective <= 0

    x_result = JuMP.value.(x; result = 1)
    @test length(x_result) == n
    @test all(isapprox(v, 0; atol = 1e-6) || isapprox(v, 1; atol = 1e-6) for v in x_result)

    raw = MOI.get(JuMP.backend(model), MOI.RawSolver())
    metadata = QUBOTools.metadata(QUBOTools.solution(raw))
    @test isempty(QUBODrivers.validate_metadata(metadata))
    @test metadata["reads"]["final_number_of_reads"] == 4
    @test metadata["seeds"]["sampler"] == 123
end

QUBODrivers.test(PySA.Optimizer; examples=true, benchmark_conformance=true) do model
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, PySA.NumberOfSweeps(), 8)
    MOI.set(model, PySA.NumberOfReplicas(), 2)
end
