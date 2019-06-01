_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/Instancias/ZRA9.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_CG_fun.jl")
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

CG_0()
#println(FPME(C))
end
# println("")
# particion()
