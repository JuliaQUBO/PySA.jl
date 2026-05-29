import PySA
import PySA: MOI, QUBODrivers
import TOML
using Test

@testset "Package metadata" begin
    project = TOML.parsefile(joinpath(pkgdir(PySA), "Project.toml"))
    qubodrivers_compat = split(project["compat"]["QUBODrivers"], r",\s*")
    readme = read(joinpath(pkgdir(PySA), "README.md"), String)

    @test "0.5" in qubodrivers_compat
    @test occursin("https://github.com/JuliaQUBO/QUBODrivers.jl", readme)
    @test !occursin("https://github.com/psrenergy/QUBODrivers.jl", readme)
end

QUBODrivers.test(PySA.Optimizer; examples=true) do model
    MOI.set(model, MOI.Silent(), true)
end
