 _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
include(_ref_*"Codigos/m_Z/parametros.jl")
include(_ref_*"Codigos/m_Z/m_Z_fun.jl")
include(_ref_*"Codigos/m_Z/m_Z_BB_fun.jl")

# #Cantidad de cuarteles
Q = sum(i for i = 1:lar) * sum(j for j = 1:anc)
#Cuarteles
C = zeros(Q, lar * anc)

#Algoritmo de generación de cuarteles
con = 0
for j = 1:anc
    for l = 0:anc-1
        if j + l <= anc
            for i = 1:lar
                for k = 0:lar-1
                    if k + i <= lar
                        global con += 1
                        for i1 = i:i+k
                            for j1 = j:j+l
                                global C[con, (i1-1)*anc+j1] = 1
                            end
                        end
                    end
                end
            end
        end
    end
end

#Matriz de varianzas
muest = muestras[:, 1]
for j = 2:anc
    global muest = [muest; muestras[:, j]]
end

zonas = Dict()
Zonas()
varianzas = zeros(Q)
for k=1:Q varianzas[k] = Vari(zonas[k]) end

BnB()

#display(Soluciones)
