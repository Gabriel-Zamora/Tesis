###############################################################################
                    # FUNCIONES #
###############################################################################
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
