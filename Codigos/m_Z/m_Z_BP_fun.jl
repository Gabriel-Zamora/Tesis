###############################################################################
                    # Branch & Price #
###############################################################################
function ZFija(Z,fijas = [0], infac = [0])
    Cand = FPM(Z,1)[:z]
    for z=1:length(Cand)-1
        if (~(Cand[z] in fijas))&(~(Cand[z] in infac))
            return Cand[z]
        end
    end
    return ""
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

function mejorcolBP_0(zfijas)
    xo = zeros(lar,anc)
    fo = 0
    col = 0
    for Col in Columnas
        dim ,H, W,cI, cJ = Col
        X = zeros(lar,anc)
        X[cI,cJ] = ones(H,W)
        f = fSP(vect,X)

        prueba = true
        for f in zfijas
            if sum(zonas[f].*X)>0 prueba = false end
        end

        if (f < -1e-5)&prueba
            if (f < fo)& probar(X)
                xo = X
                fo = f
                col = Col
            end
        end
    end
    return xo, col
end

function GenerarParticion()
    Z = [z for z=1:Q if ~(z in Infactibles)]
    zfijas = [0]
    flag = true
    Estado = :Optimo
    FO = 0

    while flag
        zfija = ZFija(Z,zfijas,Infactibles)
        if (string(termination_status(PM)) == "OPTIMAL")&(zfija!="")
            if zfijas == [0]
                zfijas = [zfija]
            else
                zfijas = [zfijas;zfija]
            end
            Z = [[z for z=Z if (sum(zonas[z].*zonas[zfija])==0)];zfija]

            global vect = FPM(Z)
            X, Col = mejorcolBP_0(zfijas)
            if sum(X) > 1
                agregarBP_0(X,Col)
                Z = [Z;Q]
            end
        else
            flag = false
        end
    end

    if string(termination_status(PM)) == "OPTIMAL"
        Estado = :Optimo
        FO = length(zfijas)
    else
        Estado = :Infactible
        FO = lar*anc
    end

    return FO,zfijas,Estado
end

function bp_0()
    flag1 = true
    while flag1
        global vect = FPM([z for z=1:Q if ~(z in Infactibles)])
        cg_0(max(lar,anc))

        ITER = iter
        flag2 = true
        while flag2
            global part = GenerarParticion()
            if part in values(Ramas)
                flag2 = false
            else
                global iter += 1
                global Ramas[iter] = part
            end
        end

        if ITER == iter
            flag1 = false
        end
    end
end

function buscaratasco(part)
    LISTA = part[2]
    optimos = Dict()
    infactibles = Dict()

    for i in LISTA
        optimos[i] = 0
        infactibles[i] = 0
    end

    for i=1:iter
        for l in LISTA
            if l in Ramas[i][2]
                if Ramas[i][3] == :Optimo
                    optimos[l] += 1
                else
                    infactibles[l] += 1
                end
            end
        end
    end

    val, key = findmin(optimos)

    if val > 0
        val2 = 0
        for i in LISTA
            if (infactibles[i] >= val2)&(optimos[i]==val)
                key = i
                val2 = infactibles[i]
            end
        end
    end

    return key
end

function BP_0(criterio = 5)
    flag = true
    cont = 1
    bp_0()
    while flag
        if cont > criterio
            flag = false
        end
        Part = part

        for i=1:cont
            global Infactibles = [Infactibles;buscaratasco(part)]
            bp_0()
        end

        global Infactibles = [0,lar*anc+1]
        bp_0()

        if part == Part
            cont +=1
        end
    end
end
