include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Z9.jl")
using Statistics, JuMP, GLPK, Gurobi

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

desabs = zeros(Q)
desi = zeros(Q)
for k=1:Q
        global desi = muest.*C[k,:]
        if sum(C[k,:])>1
            global promedio = mean(desi[i] for i=1:anc*lar if desi[i]>0)
            desabs[k] = sum(abs(desi[i]-promedio) for i=1:anc*lar if desi[i]>0)/sum(C[k,:])
        end
        global desi = zeros(Q)
end

#Parametros
a = 0.5

promedio = mean(muest)
desabst = sum(abs(muest[i]-promedio) for i=1:anc*lar)/(lar*anc)

#Modelo
# m = Model(with_optimizer(GLPK.Optimizer))
m = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0))

@variable(m, x[1:Q], Bin)

@objective(m, Min, sum(x))

@constraint(m, (1-a)*desabst*(lar*anc-sum(x)) >= sum((sum(C[i,:])-1)*desabs[i]*x[i] for i=1:Q))
@constraint(m,[j=1:lar*anc],sum(C[i,j]*x[i] for i=1:Q) == 1)

@elapsed optimize!(m)

@show objective_value(m)
@show sum([ value(x[i]) for i in 1:Q ])
