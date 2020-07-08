###############################################################################
                    # Branch & Bounds #
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

function FPR(Z=1:Q)
    Ay = Ady(Z)
    global PM = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env), "Presolve" => 0,"OutputFlag" => 0))

    @variable(PM, q[Z] >= 0)
    @variable(PM, y[Z,I,1:T] >= 0)

    @objective(PM, Min, -sum(precio[i]*rendimientos[k,i]*y[k,i,t] for t=1:T for i in I for k in Z))

    @constraint(PM,[f=1:fam-1,a in Ay,t=1:T], sum(y[k,i,t] for i in familias[f] for k in a) <= 1)
    @constraint(PM,[k=Z,t=2:T,f=1:fam-1], sum(y[k,i,τ] for i in familias[f] for τ=t-1:t) <= 1)

    @constraint(PM,rp1[k=Z,t=1:T],sum(y[k,:,t])-q[k] == 0)
    @constraint(PM,rp2[k=Z], sum(y[k,esp,:])-q[k] == 0)

    @constraint(PM,rpc,sum(q) <= L)
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z=Z) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z=Z) == 1)

    optimize!(PM)

    return objective_value(PM)
end

function FPMsE(Z=1:Q)
    Ay = Ady(Z)
    global PMsE = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env),"OutputFlag" => 0))
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
    global PMe = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env),"OutputFlag" => 0))

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

function Ady(Z)
    Ay = []
    for (i,j) in A
        if (i in Z)&(j in Z)
            Ay = [Ay; (i,j)]
        end
    end
    return Set(Ay)
end

function bnb(padre::Int64)
    if esfactible(padre)
      FPR(Nodos[padre].Z)
      global Nodos[padre].base = MBase(padre)
   else
       global Nodos[padre].base = Any[]
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

function Avanzar(padre::Int64 = numn, metodo::Symbol = :BB)
   global NP = NP[NP.Nodo .!= padre,:]

   if metodo == :BP
       bnp(padre)
   else
       bnb(padre)
   end

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
       if metodo == :BP
           global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       else
           global Nodos[numn].vect = Any[]
       end
       global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
       global Nodos[numn].base = MBase(numn)
       push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

       global numn += 1
       global Nodos[padre].suc = [Nodos[padre].suc;numn]
       global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,-objective_value(PM),[],[])]
       global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
       global Nodos[numn].var_0 = Nodos[padre].var_0
       global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
       if metodo == :BP
           global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       else
           global Nodos[numn].vect = Any[]
       end
       global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
       global Nodos[numn].base = MBase(numn)
       if (length(Nodos[numn].var_1)<L)|(Nodos[numn].var_1 == Nodos[numn].Z)
           push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])
       end

       global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)
   end

   if (metodo == :BP)&(size(NP)[1] == 0)
       act_VO()
       global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)
   end
end

function Avanzar_0(padre::Int64 = numn, metodo::Symbol = :BB)
    global NP = NP[NP.Nodo .!= padre,:]
    if (maximum(FO[:,3]) > Nodos[padre].VO)|(~esfactible(padre))
        flag = false
    elseif Nodos[padre].suc != []
        flag = false
    else
        flag = true
    end

    while flag
       global NP = NP[NP.Nodo .!= padre,:]

       if metodo == :BP
           bnp(padre)
       else
           bnb(padre)
       end

       cand = ZFija(padre,Nodos[padre].var_1,Nodos[padre].var_0)

       if cand == ""
           flag = false
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
           if metodo == :BP
               global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
           else
               global Nodos[numn].vect = Any[]
           end
           global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
           global Nodos[numn].base = MBase(numn)
           push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

           global numn += 1
           global Nodos[padre].suc = [Nodos[padre].suc;numn]
           global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,-objective_value(PM),[],[])]
           global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
           global Nodos[numn].var_0 = Nodos[padre].var_0
           global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
           if metodo == :BP
               global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
           else
               global Nodos[numn].vect = Any[]
           end
           global Nodos[numn].VO = -FPMe(Nodos[numn].Z,Nodos[numn].var_1)
           global Nodos[numn].base = MBase(numn)
           if (length(Nodos[numn].var_1)<L)|(Nodos[numn].var_1 == Nodos[numn].Z)
               push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])
           end

           global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)

           if Nodos[numn].VO > maximum(FO[:,3])
               padre = numn
           else
               flag = false
           end
       end

       if ~esfactible(padre)
           flag = false
       end
   end
   if (metodo == :BP)&(size(NP)[1] == 0)
       act_VO()
       global NP = sort!(NP[NP.FO .> maximum(FO[:,3]), :],(:FO,:NF),rev=true)
   end
end

function branching(metodo::Symbol = :BB)
   if metodo == :BP
       global vect = FPM()
   else
       global vect = Any[]
       FPR()
   end
   global numn = 1
   global Nodos = [Nodo(numn,[],[],[],[lar*anc+1],1:Q,-objective_value(PM),[],vect)]
   global Nodos[1].base = MBase(1)
   push!(NP,[1 Nodos[1].VO length(Nodos[1].var_1)])

   flag = true
   while flag
      if size(NP) == (0,3)
         flag = false
      else
          if (metodo == :BB)&(length(Soluciones) <= 2)
              Avanzar_0(NP[:,:Nodo][1],metodo)
          else
              Avanzar(NP[:,:Nodo][1],metodo)
          end
      end
   end
end

function BnB(metodo::Symbol = :BB)

   if metodo == :CG
       global Soluciones = [(0,0,FPMsE(Z),Z)]
       global soluciones = [Set(Z)]
       global FO = [0 0 objective_value(PMsE)]
       global arbol = 1
   else
       global Soluciones = [(0,0,0,Z)]
       global soluciones = [Set(Z)]
       global FO = [0 0 0]
       global arbol = 1
   end

  global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
  branching()
end
