using Statistics, JuMP, Gurobi, Plotly
###############################################################################
                    # FUNCIONES Z#
###############################################################################
function FPMZ(Z=1:Q ,modo=0)
    global PMZ = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0,gurobi_env))
    @variable(PMZ, q[Z] >= 0)
    @objective(PMZ, Min, sum(q))
    @constraint(PMZ,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z in Z) <= vt*(lar*anc)*(1-a))
    @constraint(PMZ,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z in Z) == 1)

    optimize!(PMZ)

    if modo == 0
        pii = dual(rpii)
        vpp = [dual(rpp[s]) for s=1:lar*anc]
        pp = zeros(lar,anc)
        for j=1:anc
            pp[1:lar,j] = vpp[1+lar*(j-1):lar+lar*(j-1)]
        end
        return pii,pp
    else
        zz = [0]
        vv = [0]
        su = [0]
        for z in Z
            if (value(q[z])>0)
                zz = [zz;z]
                vv = [vv;value(q[z])]
                su = [su;sum(C[z,:])]
            end
        end
        return sort(DataFrame(z=zz,val=vv,sum=su),(:val,:sum),rev=true)
    end
end


function FPMEZ(Col=C)
    Q, = size(Col)
    global PMEZ = Model(with_optimizer(Gurobi.Optimizer,Presolve=0,OutputFlag=0,gurobi_env))
    @variable(PMEZ, q[1:Q], Bin)
    @objective(PMEZ, Min, sum(q))
    @constraint(PMEZ,rpii, sum(((sum(Col[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=1:Q) <= vt*(lar*anc)*(1-a))
    @constraint(PMEZ,rpp[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PMEZ)

    return objective_value(PMEZ)
end

function fSPZ(vectz,X)
    pii = vectz[1]
    pp = vectz[2]
    return 1-sum(pp[i,j]*X[i,j] for i=1:lar for j=1:anc)-pii*((sum(X)-1)*Vari(X)+(1-a)*vt)
end


###############################################################################
                    # FUNCIONES ZRA#
###############################################################################
function FPM(Z=1:Q,Ay=A,modo=0)
    global PM = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0,gurobi_env))

    @variable(PM, q[Z] >= 0)
    @variable(PM, y[Z,I,1:T] >= 0)

    @objective(PM, Min, -sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in Z))

    @constraint(PM,[f=1:fam-1,a in Ay,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PM,[k=Z,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    global rr1 = @constraint(PM,rp1[k=Z,t=1:T],sum(y[k,:,t])-q[k] == 0)
    global rr2 = @constraint(PM,rp2[k=Z], sum(y[k,esp,:])-q[k] == 0)

    @constraint(PM,rpc,sum(q) <= L)
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=Z) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=Z) == 1)

    optimize!(PM)

    if modo==0
        pii = dual(rpii)
        pc = dual(rpc)
        vpp = [dual(rpp[s]) for s=1:lar*anc]
        pp = zeros(lar,anc)
        for j=1:anc
            pp[1:lar,j] = vpp[1+lar*(j-1):lar+lar*(j-1)]
        end

        p1 = zeros(T)
        p2 = 0
        for i=1:lar*anc
            if i in Z #hay que ver como solucionarlo
                if value(q[i]) == 0
                    p1[:] = [dual(rp1[i,t]) for t=1:T]
                    p2 = dual(rp2[i])
                end
            end
        end

        return pii,pp,pc,p1,p2
    else
        zz = []
        vv = []
        su = []
        for z in Z
            if (value(q[z])>0)
                zz = [zz;z]
                vv = [vv;value(q[z])]
                su = [su;sum(C[z,:])]
            end
        end
        return sort(DataFrame(z=zz,val=vv,sum=su),(:val,:sum),rev=true)
    end
end

function FPME(Col=C)
    Q, = size(Col)
    global PME = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))
    global q_pe = @variable(PME, q[1:Q], Bin)
    global y_pe = @variable(PME, y[1:Q,I,1:T], Bin)

    @objective(PME, Max, sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in 1:Q))

    @constraint(PME,[f=1:fam-1,a in A,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PME,[k=1:Q,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    @constraint(PME,[k=1:Q,t=1:T],sum(y[k,:,t])-q[k] == 0)
    @constraint(PME,[k=1:Q], sum(y[k,esp,:])-q[k] == 0)

    @constraint(PME,sum(q) <= L)
    @constraint(PME,sum(((sum(Col[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=1:Q) <= vt*(lar*anc)*(1-a))
    @constraint(PME,[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)

    optimize!(PME)

    return round(objective_value(PME))
end

function FPMsE(Z=1:Q)
    Ay = Ady(Z)
    global PME = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))
    global q_pe = @variable(PME, q[Z] == 1)
    global y_pe = @variable(PME, y[Z,I,1:T], Bin)

    @objective(PME, Max, sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in Z))

    @constraint(PME,[f=1:fam-1,a in Ay,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PME,[k=Z,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    @constraint(PME,[k=Z,t=1:T],sum(y[k,:,t])-q[k] == 0)
    @constraint(PME,[k=Z], sum(y[k,esp,:])-q[k] == 0)

    optimize!(PME)

    return round(objective_value(PME))
end

function fSP(vect,X)
    pii = vect[1]
    pp = vect[2]
    pc = vect[3]
    p1 = vect[4]
    p2 = vect[5]
    return sum(X)*(sum(p1)+p2)-sum(pp[i,j]*X[i,j]
    for i=1:lar for j=1:anc)-pii*((sum(X)-1)*Vari(X)+(1-a)*vt)-pc
end

function Vari(X)
    vari = muestras.*X
    if sum(X)>1
        return var(vari[i,j] for i=1:lar for j=1:anc if vari[i,j]>0 corrected=false)
    else
        return 0
    end
end

function Rendimientos(prueba = 0)
    if prueba == 1
        global rendimientos = [rendimientos;[sum(C[Q,:])*rendimiento[i] for i=1:esp]']
    else
        global rendimientos = zeros(Q,esp)
        for k=1:Q
            for i=1:esp
                global rendimientos[k,i] = sum(C[k,:])*rendimiento[i]
            end
        end
    end
end

function Adyacencia(prueba = 0)
    if prueba == 1
        for k=1:Q-1
            Ady(k,Q)
        end
    else
        for k=1:Q for l=k:Q
            Ady(k,l)
        end end
    end
end

function Zonas(prueba = 0)
    if prueba == 1
        zona = zeros(lar,anc)
        for i=1:lar
            zona[i,:] = C[Q,(i-1)*anc+1:i*anc]
        end
        global zonas[Q] = zona
    else
        for k=1:Q
            zona = zeros(lar,anc)
            for i=1:lar
                zona[i,:] = C[k,(i-1)*anc+1:i*anc]
            end
            global zonas[k] = zona
        end
    end
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

function particion()
    for t=1:T
        Matriz = zeros(lar,anc)
        for k=1:Q
            if value(q_pe[k])>0
                for i in I
                    if value(y_pe[k,i,t])>0
                        for f=1:fam
                            if i in familias[f]
                                Matriz = Matriz + f*zonas[k]
                            end
                        end
                    end
                end
            end
        end
        display(floor.(Int,Matriz))
    end
end

function mapear(t)  #Falta agregar titulos
   Matriz = zeros(lar,anc)
   iter = 0
   areas = Dict()
   for k=1:Q
       if value(q_pe[k])>0
           for i in I
               if value(y_pe[k,i,t])>0
                   for f=1:fam
                       if i in familias[f]
                           Matriz = Matriz + (f-1)*zonas[k]
                       end
                   end
               end
           end
       end
   end

   mapa = heatmap(
       x=["Z"*string(i) for i=1:lar],
       y=["Z"*string(i) for i=1:anc],
       z=Matriz,
       title = "Hola")

   layout = Layout(;title = "Periodo "*string(t))

   plot(mapa, layout)
end

###############################################################################
                    # GENERACION DE COLUMNAS Z#
###############################################################################
function agregarz_0(x, Col=0)
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
        global varianzas = [varianzas;Vari(x)]
        filter!(x->x≠Col,Columnas)
        Zonas(1)
        Rendimientos(1)
        Adyacencia(1)
        global vectz = FPMZ()
    end
end

function mejorcolCGZ_0()
    xo = zeros(lar,anc)
    fo = 0
    col = 0
    for Col in columnas
        dim ,H, W,cI, cJ = Col
        X = zeros(lar,anc)
        X[cI,cJ] = ones(H,W)
        f = fSPZ(vectz,X)
        if f < -1e-5
            if (f < fo)& probar(X)
                xo = X
                fo = f
                col = Col
            end
        elseif f > 0
            global columnas = filter(x->x≠Col,columnas)
        end
    end
    return xo, col
end

function cgz_0(niter = 1e6)
    global columnas = Columnas
    flag = true
    iter = 0
    while flag & (iter < niter)
        iter += 1
        q = Q
        X, Col = mejorcolCGZ_0()
        if sum(X) > 1
            agregarz_0(X, Col)
            FPMEZ(C)
        end
        if (q == Q)|(objective_value(PMEZ)<L)
            flag = false
        end
    end
end

function CGZ_0()
    flag = true
    while flag
        q = Q
        cgz_0()
        if (q == Q)|(objective_value(PMEZ)<L)
            flag = false
        end
    end
end

###############################################################################
                    # GENERACION DE COLUMNAS ZRA#
###############################################################################
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
        global vect = FPM()
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
        elseif f > 0
            global columnas = filter(x->x≠Col,columnas)
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
    global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
    if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

    global vectz = FPMZ()
    FPMEZ(C)
    CGZ_0()
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
