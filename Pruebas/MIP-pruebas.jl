_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Pruebas/funciones.jl")

Lista = [9]#16 25 36 49 64 81 100]

for l in 1:length(Lista)
    parametros(Lista[l])

    #BB
    tiempo = @elapsed BB.F()
    dt = CSV.read(_ref_*"/Pruebas/BB.csv",copycols=true)
    dt[:Z][l] = BB.Q
    dt[:OV][l] = objective_value(BB.m)/Lista[l]
    dt[:T][l] = tiempo
    dt[:V][l] = length(all_variables(BB.m))
    li =  list_of_constraint_types(BB.m)
    dt[:C][l] = sum(num_constraints(BB.m,li[i][1],li[i][2]) for i=1:length(li))
    CSV.write(_ref_*"/Pruebas/BB.csv",dt)

    #CG
    tiempo = @elapsed CG.F()
    dt = CSV.read(_ref_*"/Pruebas/CG.csv",copycols=true)
    dt[:Z][l] = CG.Q
    dt[:OV][l] = objective_value(CG.PME)/Lista[l]
    dt[:T][l] = tiempo
    dt[:V][l] = length(all_variables(CG.PM))
    li =  list_of_constraint_types(CG.PM)
    dt[:C][l] = sum(num_constraints(CG.PM,li[i][1],li[i][2]) for i=1:length(li))
    CSV.write(_ref_*"/Pruebas/CG.csv",dt)

    #BP
    tiempo = @elapsed BP.F()
    BP.FPM()
    dt = CSV.read(_ref_*"/Pruebas/BP.csv",copycols=true)
    dt[:Z][l] = BP.Q
    dt[:OV][l] = maximum(BP.FO[:,3])/Lista[l]
    dt[:T][l] = tiempo
    dt[:V][l] = length(all_variables(BP.PM))
    li =  list_of_constraint_types(BP.PM)
    dt[:C][l] = sum(num_constraints(BP.PM,li[i][1],li[i][2]) for i=1:length(li))
    CSV.write(_ref_*"/Pruebas/BP.csv",dt)
end
