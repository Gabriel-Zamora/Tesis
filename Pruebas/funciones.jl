using CSV, Statistics, JuMP, Gurobi, Plotly

function parametros(n)
    param = open(pwd()*"/Codigos/m_ZRA/parametros.jl")
        lineas = readlines(param)
    close(param)

    lineas[1] = "include(pwd()*\"/Codigos/Instancias/ZRA$n.jl\")"

    param = open(pwd()*"/Codigos/m_ZRA/parametros.jl","w")
        for i=1:length(lineas)
            write(param,lineas[i]*"\n")
        end
    close(param)
end

module BB
function F()
    include(pwd()*"/Codigos/m_ZRA/m_ZRA_directo.jl")
end
end

module CG
function F()
    include(pwd()*"/Codigos/m_ZRA/m_ZRA_CG.jl")
end
end
