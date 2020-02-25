push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot


N = 300
Δs = 2π/N
M = 2π

epsilon = 2e-2
potential_range = 0
rho_max = 3*M/2π

center_rho = false

atol = 1e-3
rtol = 1e-13
max_iter = 10000
step_size = 1e-2

beta_0 = 1
beta_rho0  = M/2π
beta_m     = 1
beta_h     = 1
beta_k     = 0
beta_j     = 3

if beta_h > 0
    eps_c = sqrt((beta_m^2 - beta_h/2)/beta_j^2)
else
    eps_c2 = sqrt(-beta_h/2 /beta_j^2)
end

print("Critical epsilon:")
println(eps_c)

# beta, beta_prime, beta_second = quadratic_beta(beta0, rho0, m, h)
P_3 = Bend.Params(N, Δs,
                  M, epsilon, rho_max,
                  beta_0, beta_rho0, beta_m, beta_h, beta_k, beta_j,
                  potential_range, center_rho)

P = P_3

solver_params = Bend.SolverParams(
                                  atol, rtol, max_iter, step_size,
                                  true, # adapt
                                  1e-5, # min step size
                                  1e-1, # max step size
                                  5e-1, # step down threshold
                                  5e-3, # step up threshold
                                  1.2,  # step factor
                                 )


# Circle
# Xcircle = Bend.initial_data(P, 1, 1)
# Plotting.plot(P, Xcircle, label="circle")

# Ellipsis
# Xinit = Bend.initial_data(P, 1, 1, pulse=1, poly=false, reverse_phase=false)
Xinit = Bend.initial_data_smooth(P, sides=3, smoothing=0.8)
# Xinit = Bend.initial_data(P, 1, 1, 8, true)

# Plotting.plot(P, Xinit)
# Plotting.s()

# res = Bend.minimize_energy(P, Xinit; atol=atol, max_iter=max_iter, relaxation=rel, adapt=true, dirichlet_rho=false)

# print("Tolerance met: ")
# println(res.converged)

# Plotting.plot_result(P, res, label="solution")
# Plotting.s()

Xx = copy(Xinit)

epsilons = [2e-1; 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2]
epsilons = [2e-2]

for e in epsilons
    if e == epsilons[1]
        Xx .= Xinit
    end

    f = Plotting.init()
    Plotting.plot(f, P, Xx)

    P.epsilon = e
    minimizor = Bend.minimizor(P, Xx, solver_params)

    while true
        res = minimizor()

        Xx .= res.sol
        n = res.iter
        if n % 100 == 0
            print(n)
            print("")
            print(", energy: ")
            print(res.energy_i[n])
            print(", residual norm: ")
            println(res.residual_norm_i[n])
            Plotting.plot(f, P, res.sol)
        end

        if res.finished
            print("Finished, converged: ")
            println(res.converged)
            break
        end

        # println("waiting for input")
        # readline()
    end
end

# display(X)
# println()
# println(err)
