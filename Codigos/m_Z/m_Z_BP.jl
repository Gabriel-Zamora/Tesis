 _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/Instancias/Z9.jl")
include(_ref_*"Codigos/m_Z/m_Z_CG_fun.jl")
include(_ref_*"Codigos/m_Z/m_Z_BP_fun.jl")
gurobi_env = Gurobi.Env();

#Parámetros
a = 0.5
vt = var(muestras[i,j] for i=1:lar for j=1:anc if muestras[i,j]>0 corrected=false)

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(anc*lar)']

Q, =size(C)

varianzas = zeros(lar*anc)
varianzas = [varianzas; vt]

zonas = Dict()
Zonas()

Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

Infactibles = [0,lar*anc+1]
Ramas = Dict()
iter = 0
Ramas[iter] = ([0],:Infactible,lar*anc)

# bp_0()
BP_0(2)


display(Ramas[iter])

#CG
# 102  102  129  129  129  129  129
# 102  102  129  129  129  129  129
#  63   63  195   18   96   96   96
#  63   63  195   25   96   96   96
#  63   63  195   32  146  146   35
#  63   63  195  153  153  153  153

#Directo
# 2  254  254  254  254  254  254
# 2  254  254  254  254  254  254
# 75   75   75   75  496  496  580
# 16  163  310  310  496  496  580
# 104  104  104  104  104  566  566
# 104  104  104  104  104  566  566
