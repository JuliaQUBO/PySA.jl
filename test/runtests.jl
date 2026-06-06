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

    @test "0.5" in qubodrivers_compat
    @test occursin("https://github.com/JuliaQUBO/QUBODrivers.jl", readme)
    @test !occursin("https://github.com/psrenergy/QUBODrivers.jl", readme)
    @test condapkg["pip"]["deps"]["pysa"] == "@ git+https://github.com/nasa/pysa@v0.1.0"
    @test occursin(r"version\s*=\s*v\"0\.1\.0\"", source)
end

QUBODrivers.test(PySA.Optimizer; examples=true) do model
    MOI.set(model, MOI.Silent(), true)
end
