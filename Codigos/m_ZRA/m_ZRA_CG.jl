_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZRA/parametros.jl")
include(_ref_*"Codigos/m_ZRA/m_Z_CG.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_CG_fun.jl")

C = Z_CG.C
Q, =size(C)
varianzas = Z_CG.varianzas

zonas = Dict()
Zonas()
Rendimientos()
Adyacencia()

CG_0()
println(FPME(C))


# println("")
# particion()
#
# mapear(1)
# mapear(2)
# mapear(3)
# mapear(4)
