import GLMakie
const M = GLMakie

@show M.Observable

function init_plot(P::Params, IP::IntermediateParams, S::Stiffness, X::Vector{Float64}; label::String="")

    matrices = assemble_fd_matrices(P, IP)

    fig = M.Figure(resolution = (1600, 900))

    axes = [
            M.Axis(fig[1, 1], aspect=M.AxisAspect(1)) M.Axis(fig[1, 2], title=M.L"\rho");
            M.Axis(fig[2, 1], title=M.L"\dot{\theta}, \ddot{\theta}") M.Axis(fig[2, 2], title=M.L"\beta(\rho)")
    ]

    M.limits!(axes[1,1], -1.5, 1.5, -1.5, 1.5)

    N = P.N
    Δs = IP.Δs

    t = collect(range(0, 2π, length=N+1))[1:N]

    xy = zeros(N+1, 2)

    function compute_c(_X)
        return X2candidate(P, _X; copy=false)
    end

    function compute_xy(c_node)
        for i in 2:N+1
            xy[i,:] = xy[i-1,:] + Δs*[cos(c_node.θ[i-1]); sin(c_node.θ[i-1])]
        end

        bc = sum(xy, dims=1)/N
        xy = xy .- bc

        return xy
    end

    X_node = M.Observable(X)
    c_node = M.@lift X2candidate(P, $X_node; copy=false)
    xy_node = M.@lift compute_xy($c_node)
    θ_dot_node = M.@lift compute_centered_fd_θ(P, matrices, $c_node.θ)

    # θ_dot_updown_node = M.@lift compute_fd_θ(P, matrices, $c_node.θ)

    function compute_θ_dotdot(θ)
      r = matrices.D2*θ
      r[1] -= 2π/(IP.Δs^2)
      r[end] += 2π/(IP.Δs^2)
      return r
    end

    θ_dotdot_node = M.@lift compute_θ_dotdot($c_node.θ)

    function update(X, title)
        X_node[] = X
        axes[1,1].title = title
        map(M.autolimits!, axes[2:end])
    end

    M.lines!(axes[1,1], cos.(t), sin.(t), lw=0.5, color="gray")
    M.scatterlines!(axes[1,1], M.lift(x -> x[:,1], xy_node), M.lift(x -> x[:,2], xy_node), markersize=M.lift(x -> 10*abs.([x.ρ; x.ρ[1]]), c_node), markercolor=M.lift(x -> map(v -> v>0 ? "blue" : "red", [x.ρ; x.ρ[1]]), c_node))
    M.scatter!(axes[1,1], M.lift(x -> [x[1,1]], xy_node), M.lift(x -> [x[1,2]], xy_node))
    M.lines!(axes[1,2], t, M.lift(x -> x.ρ, c_node))

    M.scatterlines!(axes[2,1], t, θ_dot_node, markersize=3, color=M.RGBAf(0.8, 0.3, 0.1, 0.5))
    M.scatterlines!(axes[2,1], t, θ_dotdot_node, markersize=3, color=M.RGBAf(0.1, 0.3, 0.8, 0.5))

    M.lines!(axes[2,2], t, M.lift(x -> S.beta(x.ρ), c_node))

    M.display(fig)

    return fig, update

end
