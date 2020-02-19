module Plotting

import PyPlot
import Bend

function init()
    PyPlot.figure(dpi=50)
end

function plot(P::Bend.Params, X::Vector{Float64}; label::String="")
    N = P.N
    Δs = P.Δs

    t = collect(range(0, 2π, length=N+1))[1:N]

    xy = zeros(N+1, 2)
    xy[2,:] = [Δs 0]

    c = Bend.X2candidate(P, X; copy=true)

    for i in 3:N+1
        xy[i,:] = xy[i-1,:] + Δs*[cos(c.θ[i-2]); sin(c.θ[i-2])]
    end

    bc = sum(xy, dims=1)/N
    xy = xy .- bc

    mass_scale = 5e0
    mass_norm = 2P.M/(Δs*mass_scale)
    mass_circles = []

    mass_circles = [mass_circles; PyPlot.matplotlib.patches.Circle(xy[1,:], c.ρ[1]/mass_norm, alpha=0.4, color="blue")]
    mass_circles = [mass_circles; PyPlot.matplotlib.patches.Circle(xy[2,:], c.ρ[2]/mass_norm, alpha=0.4, color="blue")]
    for i in 3:N
        mass_circles = [mass_circles; PyPlot.matplotlib.patches.Circle(xy[i,:], c.ρ[i]/mass_norm, alpha=0.4, color="blue")]
    end


    PyPlot.subplot(121)
    PyPlot.plot(xy[:,1], xy[:,2], label=label)
    PyPlot.gca().set_aspect("equal")
    circ_col = PyPlot.matplotlib.collections.PatchCollection(mass_circles)
    PyPlot.gca().add_collection(circ_col)
    PyPlot.gca().set_ylim(-1.5, 1.5)
    PyPlot.title("2D visualisation")

    PyPlot.subplot(222)
    PyPlot.plot(t, c.ρ, label=label)
    PyPlot.title("ρ")

    PyPlot.subplot(224)
    PyPlot.plot(t, [0; c.θ], label=label)
    PyPlot.title("θ")

    if length(label) > 0
        PyPlot.legend()
    end
end

function plot_result(P::Bend.Params, res::Bend.Result; label::String="")
    N = P.N
    Δs = P.Δs
    X = res.sol

    plot(P, X, label=label)

    PyPlot.figure(dpi=50)
    PyPlot.subplot(211)
    PyPlot.plot(res.energy_i[1:res.iter])
    PyPlot.title("Energy")
    PyPlot.subplot(212)
    PyPlot.semilogy(res.residual_i[1:res.iter])
    PyPlot.title("Residual")
end

function s()
    PyPlot.tight_layout()
    PyPlot.show()
end

end # module
