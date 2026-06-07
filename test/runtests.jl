import PySA
import PySA: MOI, QUBODrivers
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

    @test "0.5" in qubodrivers_compat
    @test occursin("https://github.com/JuliaQUBO/QUBODrivers.jl", readme)
    @test !occursin("https://github.com/psrenergy/QUBODrivers.jl", readme)
    @test occursin("## Solver options", readme)
    for attribute in (
        "PySA.NumberOfSweeps()",
        "PySA.NumberOfReplicas()",
        "PySA.NumberOfReads()",
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

@testset "README JuMP workflow" begin
    model = JuMP.Model(PySA.Optimizer)
    JuMP.set_attribute(model, MOI.Silent(), true)
    JuMP.set_attribute(model, PySA.NumberOfReads(), 4)
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
end

QUBODrivers.test(PySA.Optimizer; examples=true) do model
    MOI.set(model, MOI.Silent(), true)
end
