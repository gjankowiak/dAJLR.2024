module TestPlot

import PyPlot

function plot()
    println(PyPlot.matplotlib.patches.Circle([0; 0], 1.0))
end

end
