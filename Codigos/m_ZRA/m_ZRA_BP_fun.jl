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

        prueba = true
        for f in zfijas
            if sum(zonas[f].*X)>0 prueba = false end
        end

        if prueba
            f = fSP(vect,X)
        else
            f = 1
        end

        if (f < -1e-12)
            if (f < fo) & probar(X)
                xo = X
                fo = f
                col = Col
            end
        end
    end
    return xo, col
end

function act()
   for i=NP.Nodo
       global Nodos[i].Z = ConjZonas(Nodos[i].var_1,Nodos[i].var_0)
   end
end

function act_VO(N = numn)
   for i=1:N
        if (Nodos[i].suc == [])&(Nodos[i].var_1 != Nodos[i].Z)
            Nodos[i].Z = ConjZonas(Nodos[i].var_1,Nodos[i].var_0)
            if Nodos[i].VO != -FPMe(Nodos[i].Z,Nodos[i].var_1)
                Nodos[i].VO = -objective_value(PMe)

                if Nodos[i].VO >  maximum(FO[:,3])
                    global NP = NP[NP.Nodo .!= i,:]
                    push!(NP,[i Nodos[i].VO length(Nodos[i].var_1)])
                end

            end
        end
    end
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
            flag = false
            global Nodos[padre].Z = ConjZonas(Nodos[padre].var_1,Nodos[padre].var_0)
            if Z == Nodos[numn].Z
                act()
            end
        end
    end
end

function BnP()
   global Soluciones = [(0,0,FPMsE(Z),Z)]
   global soluciones = [Set(Z)]
   global FO = [0 0 objective_value(PMsE)]
   global arbol = 1

  global NP = DataFrame(Nodo = Int64[], FO = Float64[], NF = Int64[])
  branching(:BP)
end
