include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Datos_CEA2013.jl")
include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Pruebas/Pruebas_fun.jl")

set = [3,4,5,6,7,8,9,10]
Resultado = [:Lado :CG :CG0]

for lado=set
    println("")
    global lar = lado
    global anc = lado
    interpolar(lar,anc)
    cg = @elapsed include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Pruebas/m_Z_CG.jl")
    cg0 = @elapsed include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Pruebas/m_Z_CG_0.jl")
    println("")
    global Resultado = [Resultado; [lado cg cg0]]
    display(Resultado)
end
