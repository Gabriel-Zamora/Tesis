include(_ref_*"Codigos/Parametros/Datos_CEA2013.jl")
include(_ref_*"Codigos/Parametros/Datos_GA2018.jl")

lar = 3     #6
anc = 3     #7

L = 6

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
