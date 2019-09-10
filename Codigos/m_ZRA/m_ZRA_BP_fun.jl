###############################################################################
                    # Branch & Price #
###############################################################################
mutable struct Nodo
   ind::Int64
   pred::Array{Int64,1}
   suc::Array{Int64,1}
   var_1::Array{Int64,1}
   var_0::Array{Int64,1}
   Z::Array{Int64,1}
   VO::Float64
   base
   vect
end

function FPMsE(Z=1:Q)
    Ay = Ady(Z)
    global PMsE = Model(with_optimizer(Gurobi.Optimizer, Presolve=0,OutputFlag=0,gurobi_env))
    @variable(PMsE, y[Z,I,1:T], Bin)

    @objective(PMsE, Max, sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in Z))

    @constraint(PMsE,[f=1:fam-1,a in Ay,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PMsE,[k=Z,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    @constraint(PMsE,[k=Z,t=1:T],sum(y[k,:,t]) == 1)
    @constraint(PMsE,[k=Z], sum(y[k,esp,:]) == 1)

    optimize!(PMsE)

    return objective_value(PMsE)
end

function FPMe(Z=1:Q,QF = Z)
    Ay = Ady(Z)
    global PMe = Model(with_optimizer(Gurobi.Optimizer,OutputFlag=0,gurobi_env))

    @variable(PMe, q[Z] >= 0)
    @variable(PMe, y[Z,I,1:T] >= 0)

    @objective(PMe, Min, -sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in Z))

    @constraint(PMe,[f=1:fam-1,a in Ay,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PMe,[k=Z,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    @constraint(PMe,rp1[k=Z,t=1:T],sum(y[k,:,t])-q[k] == 0)
    @constraint(PMe,rp2[k=Z], sum(y[k,esp,:])-q[k] == 0)

    @constraint(PMe,rpc,sum(q) <= L)
    @constraint(PMe,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=Z) <= vt*(lar*anc)*(1-a))
    @constraint(PMe,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=Z) == 1)

    for k in QF for i in I for t=1:T
        delete_lower_bound(y[k,i,t])
        @constraint(PMe,y[k,i,t] in MOI.ZeroOne())
    end end end

    optimize!(PMe)
    if string(termination_status(PMe)) == "OPTIMAL"
        return objective_value(PMe)
    else
        return 0
    end
end

function ZFija(padre::Int64,fijas::Array{Int64,1} = [0], infac::Array{Int64,1} = [0])
    if Nodos[padre].base != Any[]
        Cand = Nodos[padre].base[:,1]
        for z=1:length(Cand)-1
            if (~(Cand[z] in fijas))&(~(Cand[z] in infac))
                return Cand[z]
            end
        end
    end
    return ""
end

function MBase(padre::Int64)
    Dt = DataFrame(z = Nodos[padre].Z,v = value.(all_variables(PM)[1:length(Nodos[padre].Z)]),sum = [sum(zonas[z]) for z in Nodos[padre].Z])
    Dt = sort(Dt[Dt.v .> 0, :],(:v,:sum),rev=true)
    return [Dt[:,:z] Dt[:,:v]]
end

function agregarBP_0(x, Col)
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
    end
end

function mejorcolBP_0(zfijas::Array{Int64,1},columnas = Columnas)
    xo = zeros(lar,anc)
    fo = 0
    col = 0
    for Col in columnas
        dim ,H, W,cI, cJ = Col
        X = zeros(lar,anc)
        X[cI,cJ] = ones(H,W)
        f = fSP(vect,X)

        prueba = true
        for f in zfijas
            if sum(zonas[f].*X)>0 prueba = false end
        end

        if (f < -1e-12)&prueba
            if (f < fo)& probar(X)
                xo = X
                fo = f
                col = Col
            end
        end
    end
    return xo, col
end

function ConjZonas(zfijas::Array{Int64,1},znulas::Array{Int64,1})
   prueba = [true for i=1:Q]
   for i in zfijas
       for z=1:Q
           prueba[z] &= (sum(zonas[z].*zonas[i])==0)
       end
   end

   for i in znulas
       for z=1:Q
           if i==z
           prueba[z] &= false
       end
       end
   end

   return [[z for z=1:Q if prueba[z]];zfijas]
end

function act()
   for i=1:numn
       global Nodos[i].Z = ConjZonas(Nodos[i].var_1,Nodos[i].var_0)
   end
end

function act_VO(N = numn)
   for i=1:N
        if Nodos[i].suc == []
            Nodos[i].VO = -FPMe(Nodos[i].Z,Nodos[i].var_1)
            if Nodos[i].VO >  maximum(FO[:,3])
                global NP = NP[NP.Nodo .!= i,:]
                push!(NP,[i Nodos[i].VO length(Nodos[i].var_1)])
            end
        end
    end
end

function Ady(Z)
    Ay = []
    for (i,j) in A
        if (i in Z)&(j in Z)
            Ay = [Ay; (i,j)]
        end
    end
    return Set(Ay)
end

function bnp(padre::Int64)
    if esfactible(padre)
      cg_b(padre)

      FPM(Nodos[padre].Z,Ady(Nodos[padre].Z))
      global Nodos[padre].base = MBase(padre)

   else
       global Nodos[padre].base = Any[]
   end
end

function cg_b(padre::Int64)
    flag = true
    Z = Nodos[numn].Z
    while flag
        q = Q
        global vect = FPM(Z,Ady(Z))
        X, Col = mejorcolBP_0(Nodos[padre].var_1)
        if sum(X) > 1
            agregarBP_0(X, Col)
            Z = [Z;Q]
        end
        if q == Q
            act()
            flag = false
        end
    end
end

function esfactible(nod::Int64,mod = 0)
   if mod == 0
       es = (sum(zonas[i] for i in Nodos[nod].Z) .>= ones(lar,anc))
       return ~(false in es)
   elseif mod == 1
       es = (sum(zonas[i] for i in Nodos[nod].var_1) == ones(lar,anc))
       return es
   end
end

function Avanzar(padre::Int64 = numn)
   global NP = NP[NP.Nodo .!= padre,:]

   bnp(padre)
   cand = ZFija(padre,Nodos[padre].var_1,Nodos[padre].var_0)

   if cand == ""
       if (Nodos[padre].var_1 == Nodos[padre].Z)&(~(Set(Nodos[padre].var_1) in soluciones))
           if esfactible(padre,1)&(string(termination_status(PM))== "OPTIMAL")&(maximum(FO[:,3])<FPMsE(Nodos[padre].Z))
               global Soluciones = [Soluciones;(arbol,padre,objective_value(PMsE),Nodos[padre].var_1)]
               global soluciones = [soluciones;Set(Nodos[padre].var_1)]
               global FO = [FO; [arbol padre objective_value(PMsE)]]
               global NP = sort!(NP[NP.FO .>= maximum(FO[:,3]), :],(:FO,:NF),rev=true)
           end
       end
   else
       global numn += 1
       global Nodos[padre].suc = [Nodos[padre].suc;numn]
       global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[],[cand],1:Q,-objective_value(PM),[],[])]
       global Nodos[numn].var_1 = Nodos[padre].var_1
       global Nodos[numn].var_0 = [Nodos[padre].var_0;cand]
       global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
       global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
       global Nodos[numn].base = MBase(numn)
       push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

       global numn += 1
       global Nodos[padre].suc = [Nodos[padre].suc;numn]
       global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,-objective_value(PM),[],[])]
       global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
       global Nodos[numn].var_0 = Nodos[padre].var_0
       global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
       global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
       global Nodos[numn].base = MBase(numn)
       if (length(Nodos[numn].var_1)<L)|(Nodos[numn].var_1 == Nodos[numn].Z)
           push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])
       end

       global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)
   end

   if size(NP)[1] == 0
       act_VO()
       global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)
   end
end

function branching()
   global vect = FPM()
   global numn = 1
   global Nodos = [Nodo(numn,[],[],[],[lar*anc+1],1:Q,-objective_value(PM),[],vect)]
   global Nodos[1].base = MBase(1)
   push!(NP,[1 Nodos[1].VO length(Nodos[1].var_1)])

   flag = true
   while flag
      if (size(NP) == (0,3))|(numn > 10000)
         flag = false
      else
          Avanzar(NP[:,:Nodo][1])
      end
   end
end

function BnP()
   global Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Soluciones = [(0,0,FPMsE(Z),Z)]
   global soluciones = [Set(Z)]
   global FO = [0 0 objective_value(PMsE)]
   global arbol = 1

  global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
  branching()
end
