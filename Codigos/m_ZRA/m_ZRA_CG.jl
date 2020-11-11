include(pwd()*"/Codigos/m_ZRA/parametros.jl")
include(pwd()*"/Codigos/m_ZRA/m_Z_CG.jl")
include(pwd()*"/Codigos/m_ZRA/m_ZRA_fun.jl")
include(pwd()*"/Codigos/m_ZRA/m_ZRA_CG_fun.jl")

C = Z_CG.C
Q, =size(C)
varianzas = Z_CG.varianzas
Dimensiones = Z_CG.Dimensiones
Columnas = Z_CG.Columnas

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

CG_0()

FPME(C)
#println(FPME(C))
