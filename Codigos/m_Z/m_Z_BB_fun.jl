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
    global PM = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(gurobi_env),"OutputFlag" => 0))
    @variable(PM, q[Z] >= 0)
    @objective(PM, Min, sum(q))
    @constraint(PM,rpii, sum(((sum(C[z,:])-1)*varianzas[z]+(1-a)*vt)*q[z] for z in Z) <= vt*(lar*anc)*(1-a))
    @constraint(PM,rpp[s=1:lar*anc],sum(C[z,s]*q[z] for z in Z) == 1)

    optimize!(PM)

    return objective_value(PM)
end

function ZFija(padre::Int64,fijas::Array{Int64,1} = [0], infac::Array{Int64,1} = [0])
    if Nodos[padre].base != Any[]
        Cand = Nodos[padre].base[:,1]
        for z=1:length(Cand)
            if (~(Cand[z] in fijas))&(~(Cand[z] in infac))
                return Cand[z]
            end
        end
    end
    return ""
end

function MBase(padre::Int64)
    Dt   = DataFrame(z = Nodos[padre].Z,v = value.(all_variables(PM)),sum = [sum(zonas[z]) for z in Nodos[padre].Z])
#    Dt.o = 0.5*ones(length(Dt.v))-abs.(0.5*ones(length(Dt.v))-Dt.v)
    Dt   = sort(Dt[Dt.v .> 0, :],(:v,:sum),rev=true)
    return [Dt[:,:z] Dt[:,:v]]
end

function particionBP(Lista)
   Matriz = zeros(lar,anc)
   for i in Lista
      Matriz = Matriz + i*zonas[i]
   end
   return floor.(Int,Matriz)
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
   if (minimum(FO[:,3]) < Nodos[padre].VO)|(~esfactible(padre))
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
               if esfactible(padre,1)&(string(termination_status(PM))== "OPTIMAL")
                   global Soluciones = [Soluciones;(arbol,padre,length(Nodos[padre].var_1),Nodos[padre].var_1)]
                   global soluciones = [soluciones;Set(Nodos[padre].var_1)]
                   global FO = [FO; [arbol padre length(Nodos[padre].var_1)]]
                   global NP = sort!(NP[NP.FO .<= minimum(FO[:,3]), :],(:FO,order(:NF,rev=true)))
               end
           end
       else
           global numn += 1
           global Nodos[padre].suc = [Nodos[padre].suc;numn]
           global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[],[cand],1:Q,objective_value(PM),[],[])]
           global Nodos[numn].var_1 = Nodos[padre].var_1
           global Nodos[numn].var_0 = [Nodos[padre].var_0;cand]
           global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
           if metodo == :BP
               global Nodos[numn].vect = FPM(Nodos[numn].Z)
               global Nodos[numn].VO = objective_value(PM)
           else
               global Nodos[numn].vect = Any[]
               global Nodos[numn].VO = FPR(Nodos[numn].Z)
           end
           global Nodos[numn].base = MBase(numn)
           push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

           global numn += 1
           global Nodos[padre].suc = [Nodos[padre].suc;numn]
           global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,objective_value(PM),[],[])]
           global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
           global Nodos[numn].var_0 = Nodos[padre].var_0
           global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
           if metodo == :BP
               global Nodos[numn].vect = FPM(Nodos[numn].Z)
               global Nodos[numn].VO = objective_value(PM)
           else
               global Nodos[numn].vect = Any[]
               global Nodos[numn].VO = FPR(Nodos[numn].Z)
           end
           global Nodos[numn].base = MBase(numn)
           push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

           global NP = sort!(NP[NP.FO .<= minimum(FO[:,3])-1, :],(order(:NF,rev=true),:FO))

           padre = numn
       end

       if ~esfactible(padre)
           flag = false
       end
   end

   if (metodo == :BP)&(size(NP)[1] == 0)
       act_VO()
       global NP = sort!(NP[NP.FO .<= minimum(FO[:,3])-1, :],(order(:NF,rev=true),:FO))
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
   global Nodos = [Nodo(numn,[],[],[],[lar*anc+1],1:Q,objective_value(PM),[],vect)]
   global Nodos[1].base = MBase(1)
   push!(NP,[1 Nodos[1].VO length(Nodos[1].var_1)])

   flag = true
   while flag
      if size(NP) == (0,3)
         flag = false
      else
          Avanzar(NP[:,:Nodo][1],metodo)
      end
   end
end

function BnB()
   Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Soluciones = []
   global soluciones = []
   global FO = [1 1 lar*anc+1]
   global arbol = 1

  global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
  branching()
end

#Para ZRA
function Zbranching(metodo::Symbol = :BB)
    if metodo == :BP
        global vect = FPM()
    else
        global vect = Any[]
        FPR()
    end
   global numn = 1
   global Nodos = [Nodo(numn,[],[],[],[lar*anc+1],1:Q,objective_value(PM),[],vect)]
   global Nodos[1].base = MBase(1)
   push!(NP,[1 Nodos[1].VO length(Nodos[1].var_1)])

   flag = true
   while flag
      if (size(NP) == (0,3))|(minimum(FO[:,3])<L)
         flag = false
      else
          Avanzar(NP[:,:Nodo][1])
      end
   end
end

function ZBnB()
   Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Soluciones = []
   global soluciones = []
   global FO = [1 1 lar*anc+1]
   global arbol = 0

   global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
   Zbranching()
end
