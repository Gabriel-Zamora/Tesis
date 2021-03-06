###############################################################################
                    # GENERACION DE COLUMNAS #
###############################################################################
function FPM(Z=1:Q,base=Any[])
    global PM = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "Presolve" => 0,"OutputFlag" => 0))
    @variable(PM, q[Z] >= 0)
    @objective(PM, Min, sum(q))
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z in Z) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z in Z) == 1)

    if base != Any[]
        set_start_value.(all_variables(PM), 0)

        for i in base
        set_start_value(q[i[1]],i[2])
        end
    end

    optimize!(PM)

    pii = dual(rpii)
    vpp = [dual(rpp[s]) for s=1:lar*anc]
    pp = zeros(lar,anc)
    for i=1:lar
        pp[i,1:anc] = vpp[1+anc*(i-1):anc+anc*(i-1)]
    end
    return pii,pp
end

function FPME(Col)
    Q, = size(Col)
    global PME = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "Presolve" => 0,"OutputFlag" => 0))
    global q_pe = @variable(PME, q[1:Q], Bin)
    @objective(PME, Min, sum(q))
    @constraint(PME,rpii, sum(((sum(Col[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=1:Q) <= vt*(lar*anc)*(1-a))
    @constraint(PME,rpp[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PME)

    return objective_value(PME)
end

function fSP(vect,X)
    pii = vect[1]
    pp = vect[2]
    return 1-sum(pp[i,j]*X[i,j] for i=1:lar for j=1:anc)-pii*((sum(X)-1)*Vari(X)+(1-a)*vt)
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

function agregar_0(x, Col)
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
        global vect = FPM(1:Q,[(z,value(all_variables(PM)[z])) for z in 1:Q-1])
        filter!(x->x≠Col,Columnas)
        Zonas(1)
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
    global Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
    global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
    if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]
    flag = true
    while flag
        q = Q
        cg_0()
        if q == Q
            flag = false
        end
    end
end

function Zcg_0(niter = 1e6)
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
        FPME(C)
        if (q == Q)|(objective_value(PME)<L)
            flag = false
        end
    end
end

function ZCG_0()
    global Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
    global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
    if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]
    flag = true
    while flag
        q = Q
        Zcg_0()
        if (q == Q)|(objective_value(PME)<L)
            flag = false
        end
    end
end
