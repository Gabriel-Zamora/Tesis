_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/Instancias/ZRA9.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_CG_fun.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_BP_fun.jl")
gurobi_env = Gurobi.Env();
@time begin
#Parámetros
a = 0.5
I = union(familias[1],familias[2],familias[3],familias[4])
vt = var(muestras[i,j] for i=1:lar for j=1:anc if muestras[i,j]>0 corrected=false)

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(anc*lar)']

Q, =size(C)

varianzas = zeros(lar*anc)
varianzas = [varianzas; vt]

zonas = Dict()
Zonas()
Rendimientos()
Adyacencia()

Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

vectz = FPMZ()
FPMEZ(C)
CGZ_0()

BnP(1,floor(L/2),10)

display(Soluciones)

# FPMsE(Arboles[1][51].Z)

end
#BP 9468
#CG 9814
#DD
