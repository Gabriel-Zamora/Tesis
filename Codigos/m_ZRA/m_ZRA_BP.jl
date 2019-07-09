_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZRA/parametros.jl")
include(_ref_*"Codigos/m_ZRA/m_Z_CG.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_CG_fun.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_BP_fun.jl")

C = Z_CG.C
Q, =size(C)
varianzas = Z_CG.varianzas

Z = []
for i=1:Q
    if value(all_variables(Z_CG.PME)[i]) == 1
        global Z = [Z;i]
    end
end

zonas = Dict()
Zonas()
Rendimientos()
Adyacencia()

@time begin
BnP()
end

# display(Soluciones)
# FPMsE(Arboles[1][51].Z)
