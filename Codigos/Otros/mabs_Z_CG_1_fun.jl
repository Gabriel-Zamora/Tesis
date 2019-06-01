using Statistics, JuMP, GLPK, Gurobi
gurobi_env = Gurobi.Env()

#Funciones
function FPM(Col)
    Q, = size(Col)
    global PM = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0,gurobi_env))
    @variable(PM, q[1:Q] >= 0)
    @objective(PM, Min, sum(q))
    @constraint(PM,rpii, sum(((sum(Col[z,:])-1)*desabs[z]+(1-a)*desabst)*q[z] for z=1:Q) <= desabst*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PM)

    pii = dual(rpii)
    vpp = [dual(rpp[s]) for s=1:lar*anc]
    pp = zeros(lar,anc)
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

function FSP(vect,dimen, caso = 1)
    pii = vect[1]
    pp = vect[2]
    global SP = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,  OutputFlag=0, gurobi_env))
    @variable(SP, x[1:lar,1:anc], Bin)
    @variable(SP, r[0:lar], Bin)
    @variable(SP, u[1:lar], Bin)
    @variable(SP, c[0:anc], Bin)
    @variable(SP, v[1:anc], Bin)
    @variable(SP, H >= 0, Int)
    @variable(SP, W >= 0, Int)
    global Desabs = @variable(SP, Desabs >= 0)
    @variable(SP, Promedio >= 0)
    @variable(SP, Abs[1:lar,1:anc] >= 0)

    Ez   = @expression(SP,sum(x[i,j]*muestras[i,j] for i=1:lar for j=1:anc))

    @objective(SP,Min,1-sum(pp[i,j]*x[i,j] for i=1:lar for j=1:anc)-
    pii*((dimen-1)*Desabs+(1-a)*desabst))

    @constraint(SP, Promedio <= maximum(muestras))
    @constraint(SP, Promedio*dimen == Ez)
    @constraint(SP, H*W >= dimen)

    @constraint(SP,[i=1:lar,j=1:anc],Abs[i,j] >=  muestras[i,j] - Promedio - (1-x[i,j])*maximum(muestras))
    @constraint(SP,[i=1:lar,j=1:anc],Abs[i,j] >= -muestras[i,j] + Promedio - (1-x[i,j])*maximum(muestras))
    @constraint(SP, Desabs >= sum(Abs[i,j] for i=1:lar for j=1:anc)/dimen)

    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] <= r[i])
    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] <= c[j])
    @constraint(SP,[i=1:lar,j=1:anc],x[i,j] >= r[i]+c[j]-1)
    @constraint(SP, r[0] == 0)
    @constraint(SP, c[0] == 0)
    @constraint(SP, sum(u) == 1)
    @constraint(SP, sum(v) == 1)
    @constraint(SP, sum(r) == H)
    @constraint(SP, sum(c) == W)

    if caso == 1
        @constraint(SP,[i=1:lar], r[i]-r[i-1] <= u[i])
        @constraint(SP,[j=1:anc], c[j]-c[j-1] <= v[j])
    elseif caso == 2
        @constraint(SP,[i=1:lar], r[i]-r[i-1] >= u[i])
        @constraint(SP,[j=1:anc], c[j]-c[j-1] <= v[j])
    elseif caso == 3
        @constraint(SP,[i=1:lar], r[i]-r[i-1] >= u[i])
        @constraint(SP,[j=1:anc], c[j]-c[j-1] >= v[j])
    elseif caso == 4
        @constraint(SP,[i=1:lar], r[i]-r[i-1] <= u[i])
        @constraint(SP,[j=1:anc], c[j]-c[j-1] >= v[j])
    end

    optimize!(SP)

    X = zeros(lar,anc)
    for i=1:lar for j=1:anc X[i,j] = value(x[i,j]) end end

    return objective_value(SP), X
end

function agregar(x)
    global q_s = round.(x[1,1:anc])
    for j=2:lar
        global q_s = round.([q_s; x[j,1:anc]])
    end

    global prueba = true
    for t=1:Q
        if C[t,:] == q_s global prueba = false end
    end

    if prueba
        global Q += 1
        global C = [C; q_s']
        desi = muest.*q_s
        promedio = mean(desi[i] for i=1:lar*anc if desi[i]>0)
        global desabs = [desabs ; sum(abs(desi[i]-promedio) for i=1:lar*anc if desi[i]>0)/sum(x)]
        display(round.(x))
    end
end

function inicio()
    global num = 0
    global DIMEN = Dict()
    global FO = Dict()
    global CASO = Dict()
    global X = Dict()
    global DESABS = Dict()

    vec = FPM(C)

    for dimen in Dimensiones
        for caso = 1:4
            global fo, x = FSP(vec,dimen,caso)

            global q_s = round.(x[1,1:anc])
            for j=2:lar
                global q_s = round.([q_s; x[j,1:anc]])
            end

            global prueba = true
            for t=1:Q
                if C[t,:] == q_s global prueba = false end
            end

            if sum(x) == dimen
                if prueba
                    global num += 1
                    global DIMEN[num] = dimen
                    global FO[num] = fo
                    global CASO[num] = caso
                    global X[num] = x
                    global DESABS = value(Desabs)
                end
            end
        end
    end

    for n=1:num
        if FO[n] < 0
            agregar(X[n])
        end
    end
end



function algort()
    # while
        global iter += 1

        global num = 0
        global DIMEN = Dict()
        global FO = Dict()
        global CASO = Dict()
        global X = Dict()
        global DESABS = Dict()

        vec = FPM(C)

        for dimen in Dimensiones
            for caso = 1:4
                global fo, x = FSP(vec,dimen,caso)

                global q_s = round.(x[1:lar])
                for j=2:anc
                    global q_s = round.([q_s; x[1+lar*(j-1):lar+lar*(j-1)]])
                end

                global prueba = true
                for t=1:Q
                    if C[t,:] == q_s global prueba = false end
                end

                if sum(x) == dimen
                    if prueba
                        global num += 1
                        global DIMEN[num] = dimen
                        global FO[num] = fo
                        global CASO[num] = caso
                        global X[num] = x
                        global DESABS = value(Desabs)
                    end
                end
            end
        end

        # fo, x = FSP(FPM(C),DIMEN,CASO)
        #
        # if fo<0
        #     global q_s = round.(x[1:lar])
        #     for j=2:anc
        #         global q_s = round.([q_s; x[1+lar*(j-1):lar+lar*(j-1)]])
        #     end
        #
        #     global Q += 1
        #     global C = [C; q_s']
        #     desi = muest.*q_s
        #     promedio = mean(desi[i] for i=1:lar*anc if desi[i]>0)
        #     global desabs = [desabs ; sum(abs(desi[i]-promedio) for i=1:lar*anc if desi[i]>0)]
        #     display(x)
        # else
        #     println("Terminado")
        #     break
        # end
    # end

end
