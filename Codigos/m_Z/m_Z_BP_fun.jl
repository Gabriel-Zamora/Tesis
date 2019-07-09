###############################################################################
                    # Branch & Price #
###############################################################################
mutable struct Nodo
   ind::Int16
   pred::Array{Int16,1}
   suc::Array{Int16,1}
   var_1::Array{Int16,1}
   var_0::Array{Int16,1}
   Z::Array{Int16,1}
   VO::Float64
   base
end

function ZFija(padre::Int64,fijas::Array{Int16,1} = [0], infac::Array{Int16,1} = [0])
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
    Dt = DataFrame(z = Nodos[padre].Z,v = value.(all_variables(PM)),sum = [sum(zonas[z]) for z in Nodos[padre].Z])
    Dt = sort(Dt[Dt.v .> 0, :],(:v,:sum),rev=true)
    return [Dt[:z] Dt[:v]]
end

function NBase(padre::Int64)
    base = [Nodos[padre].var_1 ones(length(Nodos[padre].var_1))]
    for i=1:length(Nodos[Nodos[padre].pred[1]].base[:,1])
        if Nodos[Nodos[padre].pred[1]].base[i,1] in setdiff(Nodos[padre].Z,Nodos[padre].var_1)
            base = [base; Nodos[Nodos[padre].pred[1]].base[i,:]']
        end
    end

    aux = zeros(lar,anc)
    for i=1:length(base[:,1])
        aux += zonas[base[i,1]]*base[i,2]
    end

    if aux != ones(lar,anc)
        for z in sort(Nodos[padre].Z,rev=true)
            if (sum(aux.*zonas[z])==0)&((sum(floor.(ones(lar,anc)-aux).*zonas[z])>0))
                aux += zonas[z]
                base = [base; [z;1]']
            end
        end
        if aux != ones(lar,anc)
            for i=1:length(base[:,1])
                if (minimum((ones(lar,anc)-aux).*zonas[base[i,1]]+(ones(lar,anc)-zonas[base[i,1]])) + base[i,2] == 1)&
                    (sum((ones(lar,anc)-aux).*zonas[base[i,1]])!=0)
                    base[i,2] = minimum((ones(lar,anc)-aux).*zonas[base[i,1]]+(ones(lar,anc)-zonas[base[i,1]]))+base[i,2]
                end
            end
            if aux != ones(lar,anc)
                for z in sort(Nodos[padre].Z,rev=true)
                    if maximum(aux.*zonas[z])<1
                        val = maximum((ones(lar,anc)-aux).*zonas[z])
                        base = [base; [z;val]']
                        aux += zonas[z]*val
                    end
                end
            end
        end
    end
    return base
end

function TBase(padre::Int64)
    base = []
    for i=1:length(Nodos[padre].base[:,1])
        base = [base; (Nodos[padre].base[i,1],Nodos[padre].base[i,2])]
    end
    return base
end

function particionBP(Lista)
   Matriz = zeros(lar,anc)
   for i in Lista
      Matriz = Matriz + i*zonas[i]
   end
   return floor.(Int,Matriz)
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
    end
end

function mejorcolBP_0(zfijas::Array{Int16,1},columnas = Columnas)
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

function ConjZonas(zfijas::Array{Int16,1},znulas::Array{Int16,1})
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

function bnp(padre::Int64)
    if esfactible(padre)
       global vect = FPM(Nodos[padre].Z,TBase(padre))
       global Nodos[padre].VO = objective_value(PM)
       global Nodos[padre].base = MBase(padre)
       X, Col = mejorcolBP_0(Nodos[padre].var_1)
       if sum(X) > 1
           agregarBP_0(X,Col)
           act()
       end
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
   if (minimum(FO[:,3]) < Nodos[padre].VO)|(~esfactible(padre))
       flag = false
   elseif Nodos[padre].suc != []
       flag = false
   else
       flag = true
   end

   while flag
       bnp(padre)

       cand = ZFija(padre,Nodos[padre].var_1,Nodos[padre].var_0)

       if (cand == "")|(minimum(FO[:,3]) < Nodos[padre].VO)
           flag = false
           if (Nodos[padre].var_1 == Nodos[padre].Z)&(~(Set(Nodos[padre].var_1) in soluciones))
               if esfactible(padre,1)&(string(termination_status(PM))== "OPTIMAL")
                   global Soluciones = [Soluciones;(arbol,padre,length(Nodos[padre].var_1),Nodos[padre].var_1)]
                   global soluciones = [soluciones;Set(Nodos[padre].var_1)]
                   global FO = [FO; [arbol padre length(Nodos[padre].var_1)]]
               end
           end
       else
           global numn += 1
           global Nodos[padre].suc = [Nodos[padre].suc;numn]
           global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[],[cand],1:Q,objective_value(PM),[])]
           global Nodos[numn].var_1 = Nodos[padre].var_1
           global Nodos[numn].var_0 = [Nodos[padre].var_0;cand]
           global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
           global Nodos[numn].base = NBase(numn)

           global numn += 1
           global Nodos[padre].suc = [Nodos[padre].suc;numn]
           global Nodos = [Nodos;Nodo(numn,[padre;Nodos[padre].pred],[],[cand],[],1:Q,objective_value(PM),[])]
           global Nodos[numn].var_1 = [cand;Nodos[padre].var_1]
           global Nodos[numn].var_0 = Nodos[padre].var_0
           global Nodos[numn].Z = ConjZonas(Nodos[numn].var_1,Nodos[numn].var_0)
           global Nodos[numn].base = NBase(numn)

           padre = numn
       end

       if ~esfactible(padre)
           flag = false
       end
   end
end

function branching()
   global vect = FPM()
   global numn = 1
   global Nodos = [Nodo(numn,[],[],[],[lar*anc+1],1:Q,objective_value(PM),[])]
   Avanzar(1)

   Lista = []
   flag = true
   while flag
      lista = [i for i=1:numn if (Nodos[i].VO <= minimum(FO[:,3]))&esfactible(i)&(Nodos[i].suc==[])&(Nodos[i].var_1 != Nodos[i].Z)]

      if (length(lista)==0)|(lista==Lista)
         flag = false
      else
         Lista = lista
         for i in lista
            Avanzar(i)
         end
      end
   end
end

function BnP()
   Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Arboles = Dict()
   global Soluciones = []
   global soluciones = []
   global FO = [1 1 lar*anc+1]
   global arbol = 0

   flag = true
   mini = lar*anc
   while flag
      global arbol += 1
      branching()
      global Arboles[arbol] = Nodos
      if mini == minimum(FO[:,3])
          FPM(1:Q)
          cg_0(lar+anc)
      end
      mini = minimum(FO[:,3])
      if arbol > 2
         if length(Arboles[arbol]) == length(Arboles[arbol-1])
            flag = false
         end
      end
   end
end
