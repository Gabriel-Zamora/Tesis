 _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/Instancias/Z9.jl")
include(_ref_*"Codigos/m_Z/m_Z_CG_fun.jl")
gurobi_env = Gurobi.Env()

#Parámetros
a = 0.5
vt = var(muestras[i,j] for i=1:lar for j=1:anc if muestras[i,j]>0 corrected=false)

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(anc*lar)']

Q, =size(C)

varianzas = zeros(lar*anc)
varianzas = [varianzas; vt]
vect = FPM(C)

zonas = Dict()
Zonas()

CG_0()
println(FPME(C))
