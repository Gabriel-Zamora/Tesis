_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZA/parametros.jl")
include(_ref_*"Codigos/m_ZA/m_Z_CG.jl")
include(_ref_*"Codigos/m_ZA/m_ZA_fun.jl")
include(_ref_*"Codigos/m_ZA/m_ZA_CG_fun.jl")
@time begin
C = Z_CG.C
Q, =size(C)
varianzas = Z_CG.varianzas

zonas = Dict()
Zonas()
Rendimientos()
Adyacencia()

CG_0()
println(FPME(C))
end
