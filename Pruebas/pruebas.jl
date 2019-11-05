_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Pruebas/funciones.jl")

Lista = [10 15 20 25 30 35 40 45 50]

for l in 6:length(Lista)
    parametros(Lista[l])

    # #BB
    # tiempo = @elapsed BB.F()
    # dt = CSV.read(_ref_*"/Pruebas/BB.csv",copycols=true)
    # dt.Z[l] = BB.Q
    # dt.OV[l] = maximum(BB.FO[:,3])/Lista[l]
    # dt.T[l] = tiempo
    # dt.V[l] = length(all_variables(BB.PM))
    # li =  list_of_constraint_types(BB.PM)
    # dt.C[l] = sum(num_constraints(BB.PM,li[i][1],li[i][2]) for i=1:length(li))
    # CSV.write(_ref_*"/Pruebas/BB.csv",dt)

    #CG
    tiempo = @elapsed CG.F()
    dt = CSV.read(_ref_*"/Pruebas/CG.csv",copycols=true)
    dt.Z[l] = CG.Q
#    dt.OV[l] = maximum(CG.FO[:,3])/Lista[l]
    dt.OV[l] = objective_value(CG.PME)/Lista[l]
    dt.T[l] = tiempo
    dt.V[l] = length(all_variables(CG.PM))
    li =  list_of_constraint_types(CG.PM)
    dt.C[l] = sum(num_constraints(CG.PM,li[i][1],li[i][2]) for i=1:length(li))
    CSV.write(_ref_*"/Pruebas/CG.csv",dt)

    #  #BP
    # tiempo = @elapsed BP.F()
    # dt = CSV.read(_ref_*"/Pruebas/BP.csv",copycols=true)
    # dt.Z[l] = BP.Q
    # dt.OV[l] = maximum(BP.FO[:,3])/Lista[l]
    # dt.T[l] = tiempo
    # dt.V[l] = length(all_variables(BP.PM))
    # li =  list_of_constraint_types(BP.PM)
    # dt.C[l] = sum(num_constraints(BP.PM,li[i][1],li[i][2]) for i=1:length(li))
    # CSV.write(_ref_*"/Pruebas/BP.csv",dt)
end

# 25,245.32975783178003,5.786520001,625,2031
# 45,245.32975783178,12.2489795,1125,5024
# 35,237.165229940191,7.2808687,875,2767
# 222,232.26651320523757,2564.172363399,5550,82365
# 210,239.88673923738733,5343.0860987,5250,63326
# 0,0,0,0,0
# 0,0,0,0,0
# 0,0,0,0,0
# 0,0,0,0,0
