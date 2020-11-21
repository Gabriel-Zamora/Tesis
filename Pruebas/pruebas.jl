include(pwd()*"/Pruebas/funciones.jl")

Lista = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100]

for l in 1:19
    parametros(Lista[l])
    # #BB
     tiempo = @elapsed BB.F()
     dt = CSV.read(pwd()*"/Pruebas/BB.csv",copycols=true)
     dt.Z[l] = BB.Q
# #    dt.OV[l] = maximum(BB.FO[:,3])/Lista[l]
     dt.OV[l] = objective_value(BB.m)/Lista[l]
     dt.T[l] = tiempo
     dt.V[l] = length(all_variables(BB.m))
     li =  list_of_constraint_types(BB.m)
     dt.C[l] = sum(num_constraints(BB.m,li[i][1],li[i][2]) for i=1:length(li))
     CSV.write(pwd()*"/Pruebas/BB.csv",dt)

    #CG
    tiempo = @elapsed CG.F()
    dt = CSV.read(pwd()*"/Pruebas/CG.csv",copycols=true)
    dt.Z[l] = CG.Q
#    dt.OV[l] = maximum(CG.FO[:,3])/Lista[l]
    dt.OV[l] = objective_value(CG.PME)/Lista[l]
#    dt.OV[l] = -objective_value(CG.PM)/Lista[l]
    dt.T[l] = tiempo
    dt.V[l] = length(all_variables(CG.PM))
    li =  list_of_constraint_types(CG.PM)
    dt.C[l] = sum(num_constraints(CG.PM,li[i][1],li[i][2]) for i=1:length(li))
    CSV.write(pwd()*"/Pruebas/CG.csv",dt)

    display(Lista[l])
    display(length(CG.A))
end
