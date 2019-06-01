include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Z9.jl")
using Statistics, JuMP, GLPK, Gurobi
gurobi_env = gurobi_env

#Parametros
a = 0.5
promedio = mean(muest)
desabst = sum(abs(muest[i]-promedio) for i=1:anc*lar)/(lar*anc)
iter = 1

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end

Q, =size(C)

muest = muestras[:,1]
for j=2:anc
    global muest = [muest; muestras[:,j]]
end

desabs = zeros(Q)

#Funciones
function FPM(C)
    global PM = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0,gurobi_env))
    @variable(PM, q[1:Q] >= 0)
    @objective(PM, Min, sum(q))
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*desabs[z]+(1-a)*desabst)*q[z] for z=1:Q) <= desabst*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PM)

    global pii = dual(rpii)
    global vpp = [dual(rpp[s]) for s=1:lar*anc]
    global pp = zeros(lar,anc)
    for j=1:anc
        global pp[1:lar,j] = vpp[1+lar*(j-1):lar+lar*(j-1)]
    end
    return pii,pp
end

function FPME(C)
    global PME = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0,gurobi_env))
    @variable(PME, q[1:Q], Bin)
    @objective(PME, Min, sum(q))
    @constraint(PME,rpii, sum(((sum(C[z,:])-1)*desabs[z]+(1-a)*desabst)*q[z] for z=1:Q) <= desabst*(lar*anc)*(1-a))
    @constraint(PME,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PME)

    return objective_value(PME)
end

function FSP(vect,H,W)
    pii = vect[1]
    pp = vect[2]
    global SP = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0,gurobi_env))
    @variable(SP, x[1:lar,1:anc], Bin)
    @variable(SP, r[0:lar], Bin)
    @variable(SP, u[1:lar], Bin)
    @variable(SP, c[0:anc], Bin)
    @variable(SP, v[1:anc], Bin)
    @variable(SP, Desabs >= 0)
    @variable(SP, Promedio >= 0)
    @variable(SP, Abs[1:lar,1:anc] >= 0)

    nz   = @expression(SP,sum(x[i,j] for i=1:lar for j=1:anc))
    Ez   = @expression(SP,sum(x[i,j]*muestras[i,j] for i=1:lar for j=1:anc))

    @objective(SP,Min,1-sum(pp[i,j]*x[i,j] for i=1:lar for j=1:anc)-
    pii*((H*W-1)*Desabs+(1-a)*desabst))

    @constraint(SP, Promedio*H*W == Ez)

    @constraint(SP,[i=1:lar,j=1:anc],Abs[i,j] >=  muestras[i,j] - Promedio - (1-x[i,j])*maximum(muestras))
    @constraint(SP,[i=1:lar,j=1:anc],Abs[i,j] >= -muestras[i,j] + Promedio - (1-x[i,j])*maximum(muestras))
    @constraint(SP, Desabs >= sum(Abs[i,j] for i=1:lar for j=1:anc))

    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] <= r[i])
    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] <= c[j])
    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] >= r[i]+c[j]-1)
    @constraint(SP, r[0] == 0)
    @constraint(SP, c[0] == 0)
    @constraint(SP, sum(u) == 1)
    @constraint(SP, sum(v) == 1)
    @constraint(SP, sum(r) == H)
    @constraint(SP, sum(c) == W)
    @constraint(SP,[i=1:lar], r[i]-r[i-1] <= u[i])
    @constraint(SP,[j=1:anc], c[j]-c[j-1] <= v[j])

    optimize!(SP)

    X = zeros(lar,anc)
    for i=1:lar for j=1:anc X[i,j] = value(x[i,j]) end end

    return objective_value(SP), round.(X)
end


    # for H=1:3
    #     for W=1:3
    #         fo, X = FSP(FPM(C),H,W)
    #
    #         global q_s = X[1:lar]
    #         for j=2:anc
    #             global q_s = [q_s; X[1+lar*(j-1):lar+lar*(j-1)]]
    #         end
    #
    #         global Q += 1
    #         global C = [C; q_s']
    #         desi = muest.*q_s
    #         promedio = mean(desi[i] for i=1:lar*anc if desi[i]>0)
    #         global desabs = [desabs ; sum(abs(desi[i]-promedio) for i=1:lar*anc if desi[i]>0)]
    #     end
    # end
