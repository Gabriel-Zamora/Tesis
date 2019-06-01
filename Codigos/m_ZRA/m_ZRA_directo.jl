include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/ZRA25.jl")
using Statistics, JuMP, GLPK, Gurobi
@time begin
#Cantidad de cuarteles
Q = sum(i for i=1:lar)*sum(j for j=1:anc)
#Cuarteles
C = zeros(Q,lar*anc)

# #Cantidad de cuarteles
Q = sum(i for i=1:lar)*sum(j for j=1:anc)
#Cuarteles
C = zeros(Q,lar*anc)

#Algoritmo de generación de cuarteles
con = 0
for j=1:anc
    for l=0:anc-1
        if j+l<= anc
            for i=1:lar
                for k=0:lar-1
                    if k+i<= lar
                        global con += 1
                        for i1=i:i+k
                            for j1=j:j+l
                                global C[con,(i1-1)*anc+j1] = 1
                            end
                        end
                    end
                end
            end
        end
    end
end

#Matriz de varianzas
muest = muestras[:,1]
for j=2:anc
    global muest = [muest; muestras[:,j]]
end

varianzas = zeros(Q)
vari = zeros(Q)
for k=1:Q
        global vari = muest.*C[k,:]
        if sum(C[k,:])>1
            varianzas[k] = var(vari[i] for i=1:anc*lar if vari[i]>0 corrected=false )
        end
        global vari = zeros(anc*lar)
end

#Matriz de rendimientos
rendimientos = zeros(Q,esp)
for k=1:Q
    for i=1:esp
        rendimientos[k,i] = sum(C[k,:])*rendimiento[i]
    end
end

#Conjunto de Adyacencias
zonas = Dict()
for k=1:Q
    zona = zeros(lar,anc)
    for i=1:lar
        zona[i,:] = C[k,(i-1)*anc+1:i*anc]
    end
    zonas[k] = zona
end

for k=1:Q for l=k:Q
    Ady(k,l)
end end


#Parametros
a = 0.5
I = union(familias[1],familias[2],familias[3],familias[4])
vt = var(muestras,corrected=false)

#Modelo
# m = Model(with_optimizer(GLPK.Optimizer))
m = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))

@variable(m, q[1:Q], Bin)
@variable(m, y[1:Q,I,1:T], Bin)

@objective(m, Max, sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in 1:Q))

@constraint(m, (1-a)*vt*(lar*anc-sum(q)) >= sum((sum(C[k,:])-1)*varianzas[k]*sum(q[k]) for k=1:Q))
@constraint(m,[j=1:lar*anc],sum(C[k,j]*sum(q[k]) for k=1:Q) == 1)
@constraint(m, sum(q) <= L)

@constraint(m,[k=1:Q,t=1:T], sum(y[k,:,t]) == q[k])

@constraint(m,[k=1:Q,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)
@constraint(m,[k=1:Q], sum(y[k,esp,:]) == q[k])

@constraint(m,[f=1:fam-1,a in A,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)

optimize!(m)

# println("")
# println("Soluciones: ")
println(objective_value(m))
# println("")
# for t=1:T
#     Matriz = zeros(lar,anc)
#     for k=1:Q
#         if value(q[k])>0
#             for i in I
#                 if value(y[k,i,t])>0
#                     for f=1:fam
#                         if i in familias[f]
#                             Matriz = Matriz + f*zonas[k]
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     display(floor.(Int,Matriz))
# end

end

#ZA
# 2x2   0.01 segundos
# 3x3   0.06 segundos
# 4x4   0.26 segundos
# 5x5   1.14 segundos
# 6x6   7.92 segundos
# 7x7  21.62 segundos
# 8x8 102.20 segundos
