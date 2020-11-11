include(pwd()*"/Codigos/Instancias/Datos_CEA2013.jl")
include(pwd()*"/Codigos/Instancias/Datos_GA2018.jl")

lar = 8     #6
anc = 8     #7

L = 13

#Interpolación
using Interpolations
itp = interpolate(muestras, BSpline(Cubic(Line(OnGrid()))))
muestras = zeros(lar,anc)
for i=0:lar-1
    for j=0:anc-1
        global muestras[i+1,j+1] = itp(1+i*5/(lar-1),1+j*6/(anc-1))
    end
end

esp = 6
fam = 4
T = 4
;
