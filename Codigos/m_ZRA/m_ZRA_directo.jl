_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZRA/parametros.jl")
include(_ref_*"Codigos/m_ZRA/m_ZRA_CG_fun.jl")

setparam!(gurobi_env, "NodefileStart", 0.5)

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

zonas = Dict()
Rendimientos()
Zonas()
Adyacencia()

#Modelo
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


# li =  list_of_constraint_types(m)
#
# display(sum(num_constraints(m,li[i][1],li[i][2]) for i=1:length(li)))
#
