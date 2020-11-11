include(pwd()*"/Codigos/Instancias/ZRA20.jl")

using Statistics, JuMP, Gurobi, Plotly
gurobi_env = Gurobi.Env();

a = 0.5
I = union(familias[1],familias[2],familias[3],familias[4])
vt = var(muestras[i,j] for i=1:lar for j=1:anc if muestras[i,j]>0 corrected=false)

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(anc*lar)']

Q, =size(C)

varianzas = zeros(lar*anc)
varianzas = [varianzas; vt]
