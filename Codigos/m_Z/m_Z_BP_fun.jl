###############################################################################
                    # Branch & Price #
###############################################################################
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

function act()
   for i=1:numn
       global Nodos[i].Z = ConjZonas(Nodos[i].var_1,Nodos[i].var_0)
   end
end

function act_VO(N = numn)
   for i=1:N
        if Nodos[i].suc == []
            FPM(Nodos[numn].Z)
            Nodos[i].VO = objective_value(PM)
            if Nodos[i].VO <=  minimum(FO[:,3])-1
                global NP = NP[NP.Nodo .!= i,:]
                push!(NP,[i Nodos[i].VO length(Nodos[i].var_1)])
            end
        end
    end
end

function bnp(padre::Int64)
    if esfactible(padre)
       cg_b(padre)

       FPM(Nodos[padre].Z)
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
        global vect = FPM(Z)
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

function BnP()
   Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Soluciones = []
   global soluciones = []
   global FO = [1 1 lar*anc+1]
   global arbol = 1

  global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
  branching(:BP)
end

#Para ZRA
function ZBnP()
   Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)
   global Columnas = [(dim,H,W,i:i+H-1,j:j+W-1) for dim in Dimensiones for H=lar:-1:1 for W=anc:-1:1
   if dim==H*W for i=1:lar-H+1 for j=1:anc-W+1]

   global Soluciones = []
   global soluciones = []
   global FO = [1 1 lar*anc+1]
   global arbol = 0

   global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
   Zbranching(:BP)
end
