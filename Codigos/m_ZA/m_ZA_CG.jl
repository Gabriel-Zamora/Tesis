_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZA/parametros.jl")
include(_ref_*"Codigos/m_ZA/m_Z_CG.jl")
include(_ref_*"Codigos/m_ZA/m_ZA_fun.jl")
include(_ref_*"Codigos/m_ZA/m_ZA_CG_fun.jl")

C = Z_CG.C
Q, =size(C)
varianzas = Z_CG.varianzas
Dimensiones = Z_CG.Dimensiones
Columnas = Z_CG.Columnas

zonas = Dict()
Zonas()
Rendimientos()
Adyacencia()

CG_0()
FPME(C)
#println(FPME(C))
