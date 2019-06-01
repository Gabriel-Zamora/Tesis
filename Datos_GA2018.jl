precio = [
0.16974
0.45
0.1923
0.42886
0.348
0
]

#Asume un cierto tamaño de suelo homogeneo
rendimiento = [
77.02962364
67.8
20
250.4544115
60.4
0
]

familias = Dict(
1 => 1:2,
2 => 3:4,
3 => 5:5,
4 => 6:6,
)


A = Set()
function Ady(n1,n2)
    z1 = zonas[n1]
    z2 = zonas[n2]
    if sum(z1.*z2) == 0
        aux = 0
        for j=1:anc-1
            aux += z1[:,j]'*z2[:,j+1]
            aux += z1[:,j+1]'*z2[:,j]
        end
        for i=1:lar-1
            aux += z1[i,:]'*z2[i+1,:]
            aux += z1[i+1,:]'*z2[i,:]
        end
        if aux > 0
            push!(A,(n1,n2))
        end
    end
end
;
