module Z_CG
_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZA/parametros.jl")
include(_ref_*"Codigos/m_Z/m_Z_CG_fun.jl")

zonas = Dict()
Zonas()

vect = FPM()
ZCG_0()
end
