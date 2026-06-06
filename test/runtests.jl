import PySA
import PySA: MOI, QUBODrivers
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
    @test pysa_pin !== nothing
    @test pysa_version !== nothing

    if pysa_pin !== nothing
        @test pysa_pin.captures[1] == "0.1.0"
    end

    if pysa_pin !== nothing && pysa_version !== nothing
        @test pysa_version.captures[1] == pysa_pin.captures[1]
    end
end

QUBODrivers.test(PySA.Optimizer; examples=true) do model
    MOI.set(model, MOI.Silent(), true)
end
