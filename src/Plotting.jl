module Plotting

import PyPlot
import Bend

function init(P::Bend.Params)
    global matrices
    matrices = Bend.assemble_fd_matrices(P)
    return PyPlot.figure(dpi=50)
end

function plot(P::Bend.Params, X::Vector{Float64}; label::String="")
    f = init(P)
    plot(f, P, X; label=label)
end

function plot(figure, P::Bend.Params, X::Vector{Float64}; label::String="")
    global matrices

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

    if length(figure.axes) > 0
        (ax1, ax2, ax3, ax4) = figure.axes
        ax1.lines[1].set_data(xy[:,1], xy[:,2])
        ax1.collections[1].set_offsets(xy)
        ax1.collections[1].set_sizes(3*c.ρ)
        ax2.lines[1].set_data(t, c.ρ)
        θ_dot = Bend.compute_centered_fd_θ(P, matrices, c.θ)
        ax3.lines[1].set_data(t, θ_dot)
        ax4.lines[1].set_data(t, Bend.compute_beta(P, c.ρ))

        ax2.relim()
        ax2.autoscale_view(true,true,true)
        ax3.relim()
        ax3.autoscale_view(true,true,true)
        ax4.relim()
        ax4.autoscale_view(true,true,true)

        PyPlot.draw()
    else
        ax1 = figure.add_subplot(221)
        ax2 = figure.add_subplot(222)
        ax3 = figure.add_subplot(223)
        ax4 = figure.add_subplot(224)

        ax1.plot(xy[:,1], xy[:,2], lw=0.1, label=label)
        ax1.scatter(xy[:,1], xy[:,2], s=3*c.ρ)
        ax1.set_aspect("equal")
        # circ_col = PyPlot.matplotlib.collections.PatchCollection(mass_circles)
        # ax1.add_collection(circ_col)
        ax1.set_xlim(-1.5, 1.5)
        ax1.set_ylim(-1.5, 1.5)
        ax1.set_title("2D visualisation")

        ax2.plot(t, c.ρ, label=label)
        ax2.axhline(0, color="black", lw=0.5)
        ax2.set_title("ρ")
        # ax2.set_ylim(0.1, 10)

        θ_dot = Bend.compute_centered_fd_θ(P, matrices, c.θ)
        ax3.plot(t, θ_dot, label=label)
        ax3.set_title("d/ds θ")
        ax3.axhline(0, color="black", lw=0.5)
        # ax3.set_ylim(0, 2π)

        ax4.plot(t, Bend.compute_beta(P, c.ρ), label=label)
        ax4.axhline(0, color="black", lw=0.5)
        ax4.set_title("β(ρ)")

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
