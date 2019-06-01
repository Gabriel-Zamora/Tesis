include(_ref_*"Codigos/Instancias/Datos_CEA2013.jl")

lar = 9     #6
anc = 9     #7

#Interpolación
using Interpolations
itp = interpolate(muestras, BSpline(Cubic(Line(OnGrid()))))
muestras = zeros(lar,anc)
for i=0:lar-1
    for j=0:anc-1
        global muestras[i+1,j+1] = itp(1+i*5/(lar-1),1+j*6/(anc-1))
    end
end

;
