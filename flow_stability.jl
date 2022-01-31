using ModulatedCurves

import Printf: @sprintf
import GLMakie

atol = 1e-9
rtol = 1e-13
max_iter = 5000
step_size = 1e-2

solver_params = SolverParams(
                                  atol, rtol, max_iter, step_size,
                                  true, # adapt
                                  1e-4, # min step size
                                  5e0, # max step size
                                  5e-1, # step down threshold
                                  1e-2, # step up threshold
                                  1.4,  # step factor
                                 )

function do_flow(P::Params, S::Stiffness,
        rho_perturbation_factor::Real, theta_smoothing::Real, perturbation_pulse::Real, do_plot::Bool=false; include_multipliers::Bool=false)

    IP = compute_intermediate(P, S)
    matrices = assemble_fd_matrices(P, IP)

    ρ_equi = P.M / P.L
    k_equi = P.L / 2π
    energy_circle = P.L * S.beta(ρ_equi) * (k_equi - P.c0)^2 / 2

    N = P.N
    t = collect(range(0, 2π, length=N + 1))[1:N]

    # perturbation = [perturbation_factor*sin.(perturbation_pulse*t).^3; perturbation_factor*sin.(perturbation_pulse*t[2:end]).^3; 0; 0; 0]

    # circ = [IP.rho_eq*ones(P.N); t[2:end]; 0; 0; 0]

    circ = ModulatedCurves.initial_data_smooth(P; sides=perturbation_pulse, smoothing=theta_smoothing)
    perturbation = [rho_perturbation_factor * sin.(perturbation_pulse * t).^3; zeros(P.N - 1); 0; 0; 0]

    Xx = circ + perturbation

    # ModulatedCurves.check_differential(P, S, Xx)
    # return
    c = ModulatedCurves.X2candidate(P, Xx)
    @show c.θ

    dθ = ModulatedCurves.compute_centered_fd_θ(P, matrices, c.θ)
    println(@sprintf "dθ = ... %.3f, %.3f, %.3f, %.3f ..." dθ[end] dθ[1] dθ[2] dθ[3])

    if do_plot
        fig, update_plot = init_plot(P, IP, S, Xx)
        # readline()
    end

    flower = build_flower(P, IP, S, Xx, solver_params; include_multipliers=include_multipliers)

    local res

    done = [false]

    iter = Iterators.takewhile((x) -> done[1], Iterators.countfrom())
    iter = 1:(max_iter-1)

    function iter_flower(_i)
        @show _i
        res = flower()

        Xx .= res.sol
        n = res.iter

        if n % 1 == 0
            println()
            print(n)

            print(", energy/rel: ")
            print(res.energy_i[n], "/", res.energy_i[n] - energy_circle)
            print(", residual norm: ")
            print(res.residual_norm_i[n])
            print(", min/max ρ: ")
            @show extrema(view(Xx, 1:P.N))
            if do_plot
                update_plot(res.sol, @sprintf "%d, energy/rel: %.4f/%.4f" n res.energy_i[n]  (res.energy_i[n] - energy_circle))
            end
        else
            print(".")
        end

        if res.finished || res.residual_norm_i[n] > 1e10
            print("Finished, converged: ")
            println(res.converged)
            @show extrema(view(Xx, 1:P.N))
            @show res.residual_norm_i[n - 1]
            @show res.residual_norm_i[n]
            if do_plot
                update_plot(res.sol, @sprintf "%d, energy/rel: %.4f/%.4f" n res.energy_i[n]  (res.energy_i[n] - energy_circle))
            end
            done[1] = true
        end
    end

    record_movie = false
    if record_movie
        GLMakie.record(iter_flower, fig, "test.mp4", iter, framerate=10)
    else
        while !done[1]
            iter_flower(1)
        end
    end
end

# constant
macro CST()
    return esc(quote
        function beta(x)
            return zero(x) .+ 1.0
        end

    function beta_prime(x)
            return zero(x)
        end

    function beta_second(x)
            return zero(x)
        end
    end)
end

# symmetrical double well
macro SDW()

    return esc(quote
        function beta(x)
            return @. 1.01 + (x - 1)^2 * (x + 1)^2
        end

    function beta_prime(x)
            return @. 4x * (x^2 - 1)
        end

    function beta_second(x)
            return @. 12x^2 - 4
        end
    end)
end

# symmetrical double well
macro SSDW()

    return esc(quote
        function beta(x)
            return @. 0.1 + (x - 1)^2 * (x + 1)^2
        end

    function beta_prime(x)
            return @. 4x * (x^2 - 1)
        end

    function beta_second(x)
            return @. 12x^2 - 4
        end
    end)
end

# asymmetrical double well
macro ADW()

    return esc(quote
        function beta(x)
            return @. 1.2 + (x - 1)^2 * (x + 1)^2 + x
        end

    function beta_prime(x)
            return @. 1 - 4x + 4x^3
        end

    function beta_second(x)
            return @. 12x^2 - 4
        end
    end)
end

function main()
    P = Params(
        # 12, # N
        3 * 4 * 5 * 3, # N
        2π,            # L
        0,             # M
        1e1,          # ε, sqrt(μ)
        1000,          # ρ_max (unused?)
        0.0,           # c0
        1, -1          # unused
        #3 * 4 * 5 * 12, 2π, 0, 1e-1, 1000, 0.0, 1, -1
    )

    # constant
    # @CST

    # symmetrical double well
    @SDW

    # asymmetrical double well
    # @ADW

    S = Stiffness(beta, beta_prime, beta_second)

    include_multipliers = false
    do_flow(P, S, 1e-1, 6e-1, 3, true; include_multipliers=include_multipliers)
end

main()
