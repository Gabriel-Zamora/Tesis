include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Z9.jl")
include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/mabs_Z_CG_1_fun.jl")

#Parametros
muest = muestras[:,1]
for j=2:anc
    global muest = [muest; muestras[:,j]]
end

a = 0.5
promedio = mean(muest)
desabst = sum(abs(muest[i]-promedio) for i=1:anc*lar)/(lar*anc)
iter = 1

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(lar*anc)']

Q, =size(C)

desabs = zeros(Q)
desabs[Q] = desabst

Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)

#Algoritmo
inicio()
