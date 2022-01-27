import GLMakie
const M = GLMakie

@show M.Observable

function init_plot(P::Params, IP::IntermediateParams, S::Stiffness, X::Vector{Float64}; label::String="")

    matrices = assemble_fd_matrices(P, IP)

    fig = M.Figure()

    axes = [
            M.Axis(fig[1, 1], aspect=M.AxisAspect(1)) M.Axis(fig[1, 2], title=M.L"\rho");
            M.Axis(fig[2, 1], title=M.L"\dot{\theta}") M.Axis(fig[2, 2], title=M.L"\beta(\rho)")
    ]

    M.limits!(axes[1,1], -1.5, 1.5, -1.5, 1.5)

    N = P.N
    Δs = IP.Δs

    t = collect(range(0, 2π, length=N+1))[1:N]

    xy = zeros(N+1, 2)
    xy[2,:] = [Δs 0]

    function compute_c(_X)
        return X2candidate(P, _X; copy=false)
    end

    function compute_xy(c_node)
        for i in 3:N+1
            xy[i,:] = xy[i-1,:] + Δs*[cos(c_node.θ[i-2]); sin(c_node.θ[i-2])]
        end

        bc = sum(xy, dims=1)/N
        xy = xy .- bc

        return xy
    end

    X_node = M.Observable(X)
    c_node = M.@lift X2candidate(P, $X_node; copy=false)
    xy_node = M.@lift compute_xy($c_node)
    θ_dot_node = M.@lift compute_centered_fd_θ(P, matrices, $c_node.θ)

    # θ_dot_min = M.@lift(

    function update(X, title)
        X_node[] = X
        axes[1,1].title = title
        map(M.autolimits!, axes[2:end])
    end

    M.lines!(axes[1,1], cos.(t), sin.(t), lw=0.5, color="gray")
    M.scatterlines!(axes[1,1], M.lift(x -> x[:,1], xy_node), M.lift(x -> x[:,2], xy_node), markersize=M.lift(x -> 10*abs.([x.ρ; x.ρ[1]]), c_node), markercolor=M.lift(x -> map(v -> v>0 ? "blue" : "red", [x.ρ; x.ρ[1]]), c_node))
    M.lines!(axes[1,2], t, M.lift(x -> x.ρ, c_node))
    M.scatterlines!(axes[2,1], t, θ_dot_node)
    M.lines!(axes[2,2], t, M.lift(x -> S.beta(x.ρ), c_node))

    M.display(fig)

    return update


    # if length(figure.axes) > 0
        # (ax1, ax2, ax3, ax4) = figure.axes
        # ax1.lines[1].set_data(xy[:,1], xy[:,2])
        # ax1.collections[1].set_offsets(xy)
        # ax1.collections[1].set_sizes(3*c.ρ)
        # ax1.collections[2].set_offsets(xy)
        # ax1.collections[2].set_sizes(-3*c.ρ)
        # ax2.lines[1].set_data(t, c.ρ)
        # θ_dot = compute_centered_fd_θ(P, matrices, c.θ)
        # ax3.lines[1].set_data(t, θ_dot)
        # ax3.lines[2].set_data(t[2:end], c.θ - t[2:end])
        # ax4.lines[1].set_data(t, S.beta(c.ρ))

        # ax2.relim()
        # ax2.autoscale_view(true,true,true)
        # ax3.relim()
        # ax3.autoscale_view(true,true,true)
        # ax4.relim()
        # ax4.autoscale_view(true,true,true)

        # PyPlot.suptitle("epsilon: " * string(P.epsilon))
        # @show P.epsilon

        # PyPlot.draw()
    # else
        # ax1 = figure.add_subplot(221)
        # ax2 = figure.add_subplot(222)
        # ax3 = figure.add_subplot(223)
        # ax4 = figure.add_subplot(224)

        # ax1.plot(xy[:,1], xy[:,2], lw=0.1, label=label)
        # ax1.scatter(xy[:,1], xy[:,2], s=3*c.ρ)
        # ax1.scatter(xy[:,1], xy[:,2], s=-3*c.ρ)
        # ax1.set_aspect("equal")
        # # circ_col = PyPlot.matplotlib.collections.PatchCollection(mass_circles)
        # # ax1.add_collection(circ_col)
        # ax1.set_xlim(-1.5, 1.5)
        # ax1.set_ylim(-1.5, 1.5)
        # ax1.set_title("2D visualisation")

        # ax2.plot(t, c.ρ, label=label)
        # ax2.axhline(0, color="black", lw=0.5)
        # ax2.set_title("ρ")
        # # ax2.set_ylim(0.1, 10)

        # θ_dot = compute_centered_fd_θ(P, matrices, c.θ)
        # ax3.plot(t, θ_dot, label=label)
        # ax3.plot(t[2:end], c.θ .- t[2:end], label=label)
        # ax3.set_title("d/ds θ")
        # ax3.axhline(0, color="black", lw=0.5)
        # # ax3.set_ylim(0, 2π)

        # ax4.plot(t, S.beta(c.ρ), label=label)
        # ax4.axhline(0, color="black", lw=0.5)
        # ax4.set_title("β(ρ)")

        # PyPlot.suptitle("epsilon: " * string(P.epsilon))

        # if length(label) > 0
            # PyPlot.legend()
        # end
end
