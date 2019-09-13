###############################################################################
                    # FUNCIONES#
###############################################################################
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
