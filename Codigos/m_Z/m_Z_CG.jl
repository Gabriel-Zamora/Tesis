include(pwd()*"/Codigos/m_Z/parametros.jl")
include(pwd()*"/Codigos/m_Z/m_Z_fun.jl")
include(pwd()*"/Codigos/m_Z/m_Z_CG_fun.jl")

zonas = Dict()
Zonas()

vect = FPM()
CG_0()

display(FPME(C))
