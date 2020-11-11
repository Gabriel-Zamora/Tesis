module Z_CG
include(pwd()*"/Codigos/m_ZRA/parametros.jl")
include(pwd()*"/Codigos/m_Z/m_Z_fun.jl")
include(pwd()*"/Codigos/m_Z/m_Z_CG_fun.jl")

zonas = Dict()
Zonas()

vect = FPM()
ZCG_0()
end
