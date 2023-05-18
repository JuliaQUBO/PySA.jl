using PySA: PySA, MOI, QUBODrivers

QUBODrivers.test(PySA.Optimizer; examples=true) do model
    MOI.set(model, MOI.Silent(), true)
end
