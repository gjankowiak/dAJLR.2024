module Plotting

import PyPlot
import Bend

function init()
    return PyPlot.figure(dpi=50)
end

function plot(P::Bend.Params, X::Vector{Float64}; label::String="")
    f = init()
    plot(f, P, X; label=label)
end

function plot(figure, P::Bend.Params, X::Vector{Float64}; label::String="")
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


    if length(figure.axes) > 0
        (ax1, ax2, ax4) = figure.axes
        ax1.lines[1].set_data(xy[:,1], xy[:,2])
        ax2.lines[1].set_data(t, c.ρ)
        ax4.lines[1].set_data(t, [0; c.θ])
        PyPlot.draw()
    else
        ax1 = figure.add_subplot(121)
        ax2 = figure.add_subplot(222)
        ax4 = figure.add_subplot(224)

        ax1.plot(xy[:,1], xy[:,2], label=label)
        ax1.set_aspect("equal")
        # circ_col = PyPlot.matplotlib.collections.PatchCollection(mass_circles)
        # ax1.add_collection(circ_col)
        ax1.set_xlim(-1.5, 1.5)
        ax1.set_ylim(-1.5, 1.5)
        ax1.set_title("2D visualisation")

        ax2.plot(t, c.ρ, label=label)
        ax2.set_title("ρ")
        ax2.set_ylim(0.1, 10)

        ax4.plot(t, [0; c.θ] .- t, label=label)
        ax4.set_title("θ")
        ax4.set_ylim(0, 2π)

        if length(label) > 0
            PyPlot.legend()
        end
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
    PyPlot.semilogy(res.residual_norm_i[1:res.iter])
    PyPlot.title("Residual norm")
end

function s()
    PyPlot.tight_layout()
    PyPlot.show()
end

end # module
