 _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_Z/parametros.jl")
include(_ref_*"Codigos/m_Z/m_Z_fun.jl")
include(_ref_*"Codigos/m_Z/m_Z_CG_fun.jl")

zonas = Dict()
Zonas()

vect = FPM()
CG_0()

display(FPME(C))
