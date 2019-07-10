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
    return [Dt[:z] Dt[:v]]
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
        global vect = Nodos[padre].vect
        X, Col = mejorcolBP_0(Nodos[padre].var_1)
       if sum(X) > 1
           agregarBP_0(X,Col)
           act()
       end

      FPM(Nodos[padre].Z,Ady(Nodos[padre].Z))
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

function Avanzar(padre::Int64 = numn)
    global NP = NP[NP.Nodo .!= padre,:]

   bnp(padre)
   cand = ZFija(padre,Nodos[padre].var_1,Nodos[padre].var_0)

   if (cand == "")|(maximum(FO[:,3])>Nodos[padre].VO)|(L<length(Nodos[padre].var_1))
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
       global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       global Nodos[numn].VO = -objective_value(PM)
       global Nodos[numn].base = MBase(numn)
       push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

       global numn += 1
       global Nodos[padre].suc = [Nodos[padre].suc;numn]
       global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,-objective_value(PM),[],[])]
       global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
       global Nodos[numn].var_0 = Nodos[padre].var_0
       global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
       global Nodos[numn].vect = FPM(Nodos[numn].Z,Ady(Nodos[numn].Z))
       global Nodos[numn].VO = -objective_value(PM)
       global Nodos[numn].base = MBase(numn)
       push!(NP,[numn Nodos[numn].VO length(Nodos[numn].var_1)])

       global NP = sort!(NP[NP.FO .>= maximum(FO[:,3]), :],(:FO,:NF),rev=true)
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
      if (size(NP) == (0,3))|(numn > sum(i for i=1:anc)*sum(i for i=1:lar)*L)
         flag = false
      else
          Avanzar(NP[:Nodo][1])
      end
   end
end

function BnP()
   global Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Arboles = Dict()
   global Soluciones = [(0,0,FPMsE(Z),Z)]
   global soluciones = [Set(Z)]
   global FO = [0 0 objective_value(PMsE)]
   global arbol = 0

   flag = true
   maxi = 0
   while flag
      global arbol += 1
      global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
      branching()
      global Arboles[arbol] = Nodos
      if maxi == maximum(FO[:,3])
          FPM(1:Q)
          cg_0(lar+anc)
      end
      maxi = maximum(FO[:,3])
      if arbol == 3 end
      if arbol > 1
         if length(Arboles[arbol])==length(Arboles[arbol-1])
            flag = false
         end
      end
   end
end
