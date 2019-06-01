using Statistics, JuMP, GLPK, Gurobi, Plotly

###############################################################################
                    # FUNCIONES #
###############################################################################
function FPM(Col)
    Q, = size(Col)
    global PM = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0,gurobi_env))
    @variable(PM, q[1:Q] >= 0)
    @objective(PM, Min, sum(q))
    @constraint(PM,rpii, sum(((sum(Col[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=1:Q) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(Col[z,s]*q[z] for z=1:Q) == 1)
    optimize!(PM)

    pii = dual(rpii)
    vpp = [dual(rpp[s]) for s=1:lar*anc]
    pp = zeros(lar,anc)
    for j=1:anc
        pp[1:lar,j] = vpp[1+lar*(j-1):lar+lar*(j-1)]
    end
    return pii,pp
end

function FPME(Col)
    Q, = size(Col)
    global PME = Model(with_optimizer(Gurobi.Optimizer,Presolve=0,OutputFlag=0,gurobi_env))
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

function Vari(X)
    vari = muestras.*X
    if sum(X)>1
        return var(vari[i,j] for i=1:lar for j=1:anc if vari[i,j]>0 corrected=false)
    else
        return 0
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
   Matriz = zeros(lar,anc)
   for i=1:Q
      if value(q_pe[i])==1
          Matriz = Matriz + i*zonas[i]
      end
   end
   return floor.(Int,Matriz)
end

function mapear()
   Matriz = zeros(lar,anc)
   iter = 0
   areas = Dict()
   for i=1:Q
      if value(q_pe[i])==1
          iter += 1
          areas[iter] = i
          Matriz = Matriz + iter*zonas[i]
      end
   end

   mapa = heatmap(
       x=["Z"*string(i) for i=1:lar],
       y=["Z"*string(i) for i=1:anc],
       z=Matriz)

   plot(mapa)

end

###############################################################################
                    # GENERACION DE COLUMNAS #
###############################################################################
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
        global vect = FPM(C)
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

###############################################################################
                    # ESTRATEGIAS 1 y 2 #
###############################################################################
function agregar_0S1(x)
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
        Zonas(1)
    end
end

function cg_0S1(niter = 1e6)
    for fm in Formas
        global vect = FPM(C)
        for Col in Columnas[fm]
            H, W, cI, cJ = Col
            X = zeros(lar,anc)
            X[cI,cJ] = ones(H,W)
            f = fSP(vect,X)
            if f < -1e-5
                agregar_0S1(X)
            end
        end
    end
end

function CG_0S1()
    Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]))
    global Columnas = Dict()
    global Formas = [(dim,H,W) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1 if dim==H*W]
    for f in Formas
        global Columnas[f] = [(f[2],f[3],i:i+f[2]-1,j:j+f[3]-1) for i=1:lar-f[2]+1 for j=1:anc-f[3]+1]
    end
    cg_0S1()
end

function CG_0S2()
    Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
    global Columnas = Dict()
    global Formas = [(dim,H,W) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1 if dim==H*W]
    for f in Formas
        global Columnas[f] = [(f[2],f[3],i:i+f[2]-1,j:j+f[3]-1) for i=1:lar-f[2]+1 for j=1:anc-f[3]+1]
    end
    cg_0S1()
end

###############################################################################
                    # ESTRATEGIA 3 #
###############################################################################
function agregar_0S3(x, Col, dim)
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
        global vect = FPM(C)
        filter!(x->x≠Col,Columnas[dim])
        Zonas(1)
    end
end

function mejorcolCG_0S3(dim)
    xo = zeros(lar,anc)
    fo = 0
    col = 0
    for Col in columnas[dim]
        H, W, cI, cJ = Col
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
            global columnas[dim] = filter(x->x≠Col,columnas[dim])
        end
    end
    return xo, col

end

function cg_0S3(niter = 1e6)
    global columnas = Columnas
    flag = true
    iter = 0
    while flag & (iter < niter)
        iter += 1
        q = Q
        for f in Formas
            X, Col = mejorcolCG_0S3(f)
            if sum(X) > 1
                agregar_0S3(X, Col, f)
            end
        end
        if q == Q
            flag = false
        end
    end
end

function CG_0S3()
    Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
    global Columnas = Dict()
    global Formas = [(dim,H,W) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1 if dim==H*W]
    for f in Formas
        global Columnas[f] = [(f[2],f[3],i:i+f[2]-1,j:j+f[3]-1) for i=1:lar-f[2]+1 for j=1:anc-f[3]+1]
    end
    flag = true
    while flag
        q = Q
        cg_0S3()
        if q == Q
            flag = false
        end
    end
end

###############################################################################
                    # ESTRATEGIA 4 #
###############################################################################
function agregar_0S4(x, Col, dim)
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
        filter!(x->x≠Col,Columnas[dim])
        Zonas(1)
    end
end

function mejorcolCG_0S4(dim)
    xo = Dict()
    fo = Dict()
    col = Dict()
    iter = 0
    for Col in columnas[dim]
        H, W, cI, cJ = Col
        X = zeros(lar,anc)
        X[cI,cJ] = ones(H,W)
        f = fSP(vect,X)
        if f < -1e-5
            if probar(X)
                iter += 1
                xo[iter] = X
                fo[iter] = f
                col[iter] = Col
            end
        elseif f > 0
            global columnas[dim] = filter(x->x≠Col,columnas[dim])
        end
    end
    return xo, fo, col

end

function cg_0S4(niter = 1e6)
    global columnas = Columnas
    flag = true
    iter = 0
    while flag & (iter < niter)
        iter += 1
        q = Q
        for f in Formas
            X, FO ,Col = mejorcolCG_0S4(f)  #Aun no defino como escoger los mejores k
            candidatos = sort([(FO[i],i) for i=1:length(FO)])
            if K < length(FO)
                for i=1:K
                    if sum(X[candidatos[i][2]]) > 1
                        agregar_0S4(X[candidatos[i][2]], Col[candidatos[i][2]], f)
                    end
                end
            else
                for i=1:length(FO)
                    if sum(X[candidatos[i][2]]) > 1
                        agregar_0S4(X[candidatos[i][2]], Col[candidatos[i][2]], f)
                    end
                end
            end
            global vect = FPM(C)
        end
        if q == Q
            flag = false
        end
    end
end

function CG_0S4()
    global K = 15
    Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
    global Columnas = Dict()
    global Formas = [(dim,H,W) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1 if dim==H*W]
    for f in Formas
        global Columnas[f] = [(f[2],f[3],i:i+f[2]-1,j:j+f[3]-1) for i=1:lar-f[2]+1 for j=1:anc-f[3]+1]
    end
    flag = true
    while flag
        q = Q
        cg_0S4()
        if q == Q
            flag = false
        end
    end
end