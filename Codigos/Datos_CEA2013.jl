using DataFrames, CSV
datos = CSV.read("C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/Datos_CEA2013.csv")

#Matriz con la muestras
muestras = zeros(6,7)
for j=1:7
    for i=1:6
        if i+(j-1)*6 <= 40
            muestras[i,j] = datos[:OM][(j-1)*6+i]
        end
    end
end
muestras[2:5,7] = muestras[1:4,7]
muestras[1,7] = 0.5*muestras[1,6]+0.5*muestras[2,7]
muestras[6,7] = 0.5*muestras[6,6]+0.5*muestras[5,7]
;
