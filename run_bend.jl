push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot
import Serialization


N = 2*2*3*5*7
Δs = 2π/N
M = 2π

epsilon = 2e-2
potential_range = 0
rho_max = 3*M/2π

center_rho = false

atol = 1e-8
rtol = 1e-13
max_iter = 50000
step_size = 1e-2

beta_0 = 1
beta_rho0  = M/2π
beta_m     = 1
beta_h     = 1
beta_k     = 0
beta_j     = 4

print("Critical epsilons:")
for j in 2:4
    if beta_h > 0
        eps_c = sqrt((beta_m^2 - beta_h/2)/j^2)
    else
        eps_c = sqrt(-beta_h/2 /j^2)
    end
    print(eps_c)
    print(" ")
end
println()

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
                                  5e1, # step down threshold
                                  5e-3, # step up threshold
                                  1.4,  # step factor
                                 )


# Circle
# Xcircle = Bend.initial_data(P, 1, 1)
# Plotting.plot(P, Xcircle, label="circle")

# Ellipsis
# Xinit = Bend.initial_data(P, 1, 1, pulse=2, reverse_phase=true)
# Xinit = Bend.initial_data_smooth(P, sides=2, smoothing=0.9, reverse_phase=true)
Xinit = Serialization.deserialize("xstart")

Xx = copy(Xinit)

epsilons = [range(0.35, 0.21, length=50); 2e-1; 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2]
# epsilons = [0.35; 0.345; 0.34]
# epsilons = [1.2455e0]
# close to critical for (m,h,k) = (-2,-2,20)
# epsilons = [1.364625e0]

# epsilons = [0.35]
# epsilons = []

Xs = []

for e in epsilons
    if e == epsilons[1]
        Xx .= Xinit
    end

    # f = Plotting.init(P)
    # Plotting.plot(f, P, Xx)

    P.epsilon = e
    minimizor = Bend.minimizor(P, Xx, solver_params)

    local res
    global Xs

    while true
        res = minimizor()

        Xx .= res.sol
        n = res.iter
        # if n % 100 == 0
            # println()
            # print(n)
            # print("")
            # print(", energy: ")
            # print(res.energy_i[n])
            # print(", residual norm: ")
            # println(res.residual_norm_i[n])
            # Plotting.plot(f, P, res.sol)
        # else
            # print(".")
        # end
        print(".")

        if res.finished
            print("Finished, converged: ")
            println(res.converged)
            push!(Xs, res.sol)
            # if e == epsilons[1]
                # Serialization.serialize("xstart", res.sol)
            # end
            break
        end

        # println("waiting for input")
        # readline()
    end

    # Plotting.plot_result(P, res, label="solution")
    # Plotting.s()
end

Serialization.serialize("branch_1.dat", Xs)
Serialization.serialize("branch_1_P.dat", P)

for r in Xs
    Plotting.plot(P, r)
end

# display(X)
# println()
# println(err)
