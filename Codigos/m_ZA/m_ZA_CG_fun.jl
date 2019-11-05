###############################################################################
                    # GENERACION DE COLUMNAS#
###############################################################################
function FPM(Z=1:Q,Ay=A,base=Any[])
    global PM = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0,gurobi_env))

    @variable(PM, q[Z] >= 0)
    @variable(PM, y[Z,I] >= 0)

    @objective(PM, Min, -sum(precio[i]*rendimientos[k,i]*y[k,i] for i in I for k in Z))

    @constraint(PM,[f=1:fam-1,a in Ay], sum(y[k,i] for i in familias[f] for k in a) <= 1)

    @constraint(PM,rp[k=Z],sum(y[k,:])-q[k] == 0)

    @constraint(PM,rpc,sum(q) <= L)
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=Z) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=Z) == 1)

    if base != Any[]
        set_start_value.(all_variables(PM), 0)

        for i in base
        set_start_value(q[i[1]],i[2])
        end
    end

    optimize!(PM)

    global cA = [[] []]
    for z in Z
        for f=1:fam-1
            if sum(value(y[z,i]) for i=familias[f]) == 1
                global cA = [cA;z f]
            end
        end
    end

    pii = dual(rpii)
    pc = dual(rpc)
    vpp = [dual(rpp[s]) for s=1:lar*anc]
    pp = zeros(lar,anc)
    for i=1:lar
        pp[i,1:anc] = vpp[1+anc*(i-1):anc+anc*(i-1)]
    end

    return pii,pp,pc,cA
end

function FPME(Col=C)
    Q, = size(Col)
    global PME = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))
    @variable(PME, q[1:Q], Bin)
    @variable(PME, y[1:Q,I], Bin)

    @objective(PME, Max, sum(precio[i]*rendimientos[k,i]*y[k,i] for i in I for k in 1:Q))

    @constraint(PME,[f=1:fam-1,a in A], sum(y[k,i] for i in familias[f] for k in a) <= 1)

    @constraint(PME,[k=1:Q],sum(y[k,:])-q[k] == 0)

    @constraint(PME,sum(q) <= L)
    @constraint(PME,sum(((sum(Col[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=1:Q) <= vt*(lar*anc)*(1-a))
    @constraint(PME,[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)

    optimize!(PME)

    return objective_value(PME)
end

function fSP(vect,X)
    pii = vect[1]
    pp = vect[2]
    pc = vect[3]
    cA = vect[4]

    Rq = - sum(pp[i,j]*X[i,j] for i=1:lar for j=1:anc) - pii*((sum(X)-1)*Vari(X)+(1-a)*vt) - pc
    return FsA(X,cA,Rq) + Rq
end

function FsA(X,cA,Rq)
    sRA = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))
    global y = @variable(sRA, y[I], Bin)

    @objective(sRA, Min, -sum(precio[i]*rendimiento[i]*sum(X)*y[i] for i in I))

    @constraint(sRA, sum(y) == 1)

    optimize!(sRA)

    return objective_value(sRA)
end

function probar(x)
    q_s =  round.(x[1,1:anc])
    for j=2:lar
        q_s = round.([q_s; x[j,1:anc]])
    end

    prueba = true
    for t=1:Q
        if C[t,:] == q_s
            prueba = false
        end
    end
    return prueba
end

function agregar_0(x, Col=0)
    q_s = round.(x[1,1:anc])
    for j=2:lar
        q_s = round.([q_s; x[j,1:anc]])
    end

    prueba = true
    for t=1:Q
        if C[t,:] == q_s prueba = false end
    end

    if prueba
        global Q += 1
        global C = [C; q_s']
        global varianzas = [varianzas ; Vari(x)]
        filter!(x->x≠Col,Columnas)
        Zonas(1)
        Rendimientos(1)
        Adyacencia(1)
        global vect = FPM(1:Q,A,[(z,value(all_variables(PM)[z])) for z in 1:Q-1])
    end
end

function mejorcolCG_0()
    xo = zeros(lar,anc)
    fo = 0
    col = 0
    for Col in columnas
        dim ,H, W,cI, cJ = Col
        X = zeros(lar,anc)
        X[cI,cJ] = ones(H,W)
        f = fSP(vect,X)
        if f < -1e-5
            if (f < fo)& probar(X)
                xo = X
                fo = f
                col = Col
            end
        end
    end
    return xo, col
end

function cg_0(niter = 1e6)
    global columnas = Columnas
    flag = true
    iter = 0
    while flag & (iter < niter)
        iter += 1
        q = Q
        X, Col = mejorcolCG_0()
        if sum(X) > 1
            agregar_0(X, Col)
        end
        if q == Q
            flag = false
        end
    end
end

function CG_0()
    global vect = FPM()

    flag = true
    while flag
        q = Q
        cg_0()
        if q == Q
            flag = false
        end
    end
end
