using CSV, Statistics, JuMP, Gurobi, Plotly

function parametros(n)
    param = open(_ref_*"Codigos/m_ZA/parametros.jl")
        lineas = readlines(param)
    close(param)

    lineas[2] = "include(_ref_*\"Codigos/Instancias/ZRA$n.jl\")"

    param = open(_ref_*"Codigos/m_ZA/parametros.jl","w")
        for i=1:length(lineas)
            write(param,lineas[i]*"\n")
        end
    close(param)
end

module BB
function F()
    _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
    include(_ref_*"Codigos/m_ZA/m_ZA_directo.jl")
end
end

module CG
function F()
    _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
    include(_ref_*"Codigos/m_ZA/m_ZA_CG.jl")
end
end

module BP
function F()
    _ref_ = "C:/Users/gz_am/Dropbox/u/Proyecto de Tésis/Julia JuMP 0.19/Tesis/"
    include(_ref_*"Codigos/m_ZA/m_ZA_BP.jl")
end
end
