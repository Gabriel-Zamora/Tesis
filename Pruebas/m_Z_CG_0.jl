include("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/m_Z_CG_0_fun.jl")
#Parámetros
a = 0.5
vt = var(muestras[i,j] for i=1:lar for j=1:anc if muestras[i,j]>0 corrected=false)

C = zeros(anc*lar,lar*anc)
for j=1:anc*lar global C[j,j] = 1 end
C = [C; ones(anc*lar)']

Q, =size(C)

varianzas = zeros(lar*anc)
varianzas = [varianzas; vt]

Dimensiones = sort(unique([i*j for i=1:max(lar,anc) for j=1:min(lar,anc) if (i*j != 1)&(i*j != lar*anc)]),rev=true)

flag = true
while flag
    q = Q

    for dim in Dimensiones
        for H=lar:-1:1 for W=anc:-1:1
            if dim == H*W
                global vect = FPM(C)
                X = mejorcol(H,W)
                if sum(X) > 1
                    agregar(X)
                end
            end
        end end
    end

    if q == Q
        global flag = false
    end
end

println("CG0 ",FPME(C))
