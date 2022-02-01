using ModulatedCurves

import Printf: @sprintf

function build_initial_data(P::Params, rho_perturbation_factor::Real, theta_smoothing::Real, perturbation_pulse::Real)
    N = P.N
    t = collect(range(0, 2π, length=N + 1))[1:N]
    circ = ModulatedCurves.initial_data_smooth(P; sides=perturbation_pulse, smoothing=theta_smoothing)
    perturbation = [rho_perturbation_factor * sin.(perturbation_pulse * t).^3; zeros(P.N); 0; 0; 0]

    return circ + perturbation
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
            return @. 1.7 + (x - 1)^2 * (x + 1)^2 + x
        end

    function beta_prime(x)
            return @. 1 - 4x + 4x^3
        end

    function beta_second(x)
            return @. 12x^2 - 4
        end
    end)
end

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

function main()
    P = Params(
        # 12, # N
        3 * 4 * 5 * 3, # N
        2π,            # L
        0,             # M
        3.3e-1,          # ε, sqrt(μ)
        1000,          # ρ_max (unused?)
        0.0,           # c0
        1, -1          # unused
        #3 * 4 * 5 * 12, 2π, 0, 1e-1, 1000, 0.0, 1, -1
    )

    # constant
    # @CST

    # symmetrical double well
    # @SDW

    # asymmetrical double well
    @ADW

    S = Stiffness(beta, beta_prime, beta_second)

    Xinit = build_initial_data(P, 1e-1, 8e-1, 3)

    include_multipliers = false
    do_flow(P, S, Xinit, solver_params; do_plot=true, include_multipliers=include_multipliers, record_movie=false)
end

main()
