include(pwd()*"/Codigos/m_Z/parametros.jl")
include(pwd()*"/Codigos/m_Z/m_Z_fun.jl")

setparam!(gurobi_env, "NodefileStart", 0.5)

# #Cantidad de cuarteles
Q = sum(i for i = 1:lar) * sum(j for j = 1:anc)
#Cuarteles
C = zeros(Q, lar * anc)

#Algoritmo de generación de cuarteles
con = 0
for j = 1:anc
    for l = 0:anc-1
        if j + l <= anc
            for i = 1:lar
                for k = 0:lar-1
                    if k + i <= lar
                        global con += 1
                        for i1 = i:i+k
                            for j1 = j:j+l
                                global C[con, (i1-1)*anc+j1] = 1
                            end
                        end
                    end
                end
            end
        end
    end
end

#Matriz de varianzas
muest = muestras[:, 1]
for j = 2:anc
    global muest = [muest; muestras[:, j]]
end

zonas = Dict()
Zonas()
varianzas = zeros(Q)
for k=1:Q varianzas[k] = Vari(zonas[k]) end

#Modelo
# m = Model(with_optimizer(GLPK.Optimizer))
m = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "Presolve" => 0,"OutputFlag" => 0))

@variable(m, q[1:Q], Bin)

@objective(m, Min, sum(q))

@constraint( m, (1 - a) * vt * (lar * anc - sum(q)) >= sum((sum(C[i, :]) - 1) *
            varianzas[i] * q[i] for i = 1:Q))
@constraint(m, [j = 1:lar*anc], sum(C[i, j] * q[i] for i = 1:Q) == 1)

optimize!(m)

display(objective_value(m))
