using PySA
using JuMP

PySA.test(; examples=true) do model
    JuMP.set_silent(model)
end
