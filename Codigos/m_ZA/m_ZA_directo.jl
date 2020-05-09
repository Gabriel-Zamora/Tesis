_ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_ZA/parametros.jl")
include(_ref_*"Codigos/m_ZA/m_ZA_fun.jl")

setparam!(gurobi_env, "NodefileStart", 0.1)
#setparam!(gurobi_env, "Threads", 1)

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

zonas = Dict()
Rendimientos()
Zonas()
Adyacencia()
varianzas = zeros(Q)
for k=1:Q varianzas[k] = Vari(zonas[k]) end

#Modelo
#m = Model(with_optimizer(Gurobi.Optimizer,Presolve=0,OutputFlag=0,gurobi_env))
m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "Presolve" => 0,"OutputFlag" => 0))
set_time_limit_sec(m, 20000)

@variable(m, q[1:Q], Bin)
@variable(m, y[1:Q,I], Bin)

@objective(m, Max, sum(precio[i]*rendimientos[k,i]*y[k,i] for i in I for k in 1:Q))

@constraint(m, (1-a)*vt*(lar*anc-sum(q)) >= sum((sum(C[k,:])-1)*varianzas[k]*sum(q[k]) for k=1:Q))
@constraint(m,[j=1:lar*anc],sum(C[k,j]*sum(q[k]) for k=1:Q) == 1)
@constraint(m, sum(q) <= L)

@constraint(m,[k=1:Q], sum(y[k,:]) == q[k])

@constraint(m,[f=1:fam-1,a in A], sum(y[k,i] for i in familias[f] for k in a) <= 1)

optimize!(m)

#display(objective_value(m))
