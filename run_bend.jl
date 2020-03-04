push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot
import Serialization


N = 2*2*3*5*7
N = 2*2*2*2*2*2*7
print("N: ")
println(N)
Δs = 2π/N
M = 2π

epsilon = 2e-2
potential_range = 0.1
rho_max = -1e-12

center_rho = false

atol = 1e-8
rtol = 1e-13
max_iter = 3000
step_size = 1e-2

beta_0 = 1
beta_rho0  = M/2π
beta_m     = -1.3
beta_h     = 1
beta_k     = 0
beta_j     = 4

print("Critical epsilons:")
for j in 2:8
    eps_c = sqrt((beta_m^2 - beta_h/2)/j^2)
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
                                  5e-4, # step up threshold
                                  1.4,  # step factor
                                 )


# Circle
# Xcircle = Bend.initial_data(P, 1, 1)
# Plotting.plot(P, Xcircle, label="circle")

# Ellipsis
# Xinit = Bend.initial_data(P, 1, 1, pulse=8, pulse_amplitude=1e-1, reverse_phase=true)
Xinit = Bend.initial_data_smooth(P, sides=8, smoothing=0.4, reverse_phase=true)
#
Xinit = Serialization.deserialize("case_6_branch_1_new.dat.bak")[1]
# res_prev = Serialization.deserialize("results_case_1/branch_2.dat")
# Xinit = res_prev[end]

Xx = Base.copy(Xinit)

# bifurcation_1
# epsilons = [range(0.35, 0.21, length=50); 2e-1; 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2]

# bifurcation_2
# epsilons = [range(0.2355, 0.20, length=60); 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2][1:2]
# epsilons = [0.2355; 0.23545; 0.23525; range(0.235, 0.2, step=-0.0005); 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2]
# epsilons = [range(0.2, 0.1, length=10); 9e-2; 8e-2; 7e-2; 6e-2; 5e-2]
# epsilons = 5*10. .^range(-2, -4; length=10)[1:6]
# println(size(epsilons))
# epsilons = [0.35; 0.345; 0.34]
# epsilons = [1.2455e0]
# close to critical for (m,h,k) = (-2,-2,20)
# epsilons = [1.364625e0]

# Xstored = Serialization.deserialize("branch_1.dat")
# Xinit .= Xstored[end]

epsilons = [1.38835]
epsilons = [1.3883]
# Case 2, oscillations
epsilons = [1.47]
# Case 2, 1.5 looks critical
# epsilons = 1 ./exp.(range(log(1/0.75), log(1/0.97), length=20))
# epsilons = [0.75; 0.755; 0.76; 0.765; 0.77; 0.775; range(0.78, 0.90, step=0.01); range(0.901, 0.95, step=0.002)]
#
# epsilons = [1.1185; 1.119; range(1.12, 1.13, step=0.001); 1 ./exp.(range(log(1/1.135), log(1/1.466), length=20))]
#
# case 3
# epsilons = [0.3]
# epsilons = 1 ./exp.(range(log(1/0.40), log(1/0.4082), length=10))
#
# case 4
epsilons = [range(0.353, 0.35, step=-0.0005); range(0.345, 0.33, step=-0.005); range(0.32, 0.05, step=-0.01)]
epsilons = collect(range(0.05, 0.01, step=-0.002))
# epsilons = [0.353]

# case 5
epsilons = [range(0.353, 0.35, step=-0.0005);
            range(0.345, 0.33, step=-0.005);
            range(0.32, 0.05, step=-0.01);
            range(0.05, 0.01, step=-0.002)]
epsilons = epsilons[59-38+1:end]
#
# case 6
epsilons = [0.136]
# epsilons = collect(range(0.136, 0.135, step=-0.0001))
# epsilons = collect(range(0.135, 0.13, step=-0.001))
# epsilons = collect(range(0.13, 0.05, step=-0.01))
epsilons = collect(range(0.05, 0.01, step=-0.005))


@show epsilons

Xs = []
Ps = []

f = Plotting.init(P)
Plotting.plot(f, P, Xx)

for e in epsilons
    println("##")
    println(e)
    println("##")

    if e == epsilons[1]
        Xx .= Xinit
    end

    P.epsilon = e
    println(P.epsilon)

    minimizor = Bend.minimizor(P, Xx, solver_params)

    local res
    global Xs

    while true
        res = minimizor()

        Xx .= res.sol
        n = res.iter
        if n % 100 == 0
            println()
            print(n)
            print("")
            print(", energy: ")
            print(res.energy_i[n])
            print(", residual norm: ")
            print(res.residual_norm_i[n])
            print(", min/max ρ: ")
            println(extrema(view(Xx, 1:P.N)))
            Plotting.plot(f, P, res.sol)
        else
            print(".")
        end
        # print(".")

        if res.finished
            print("Finished, converged: ")
            println(res.converged)
            push!(Xs, res.sol)
            push!(Ps, Bend.copy(P))
            # if e == epsilons[1]
                # Serialization.serialize("xstart_2.dat", res.sol)
            # end
            break
        end
    end

    # Plotting.plot_result(P, res, label="solution")
    # Plotting.s()
end

Serialization.serialize("case_6_branch_1.dat", Xs)
Serialization.serialize("case_6_branch_1_P.dat", Ps)

# for r in Xs
    # Plotting.plot(P, r)
# end

# display(X)
# println()
# println(err)
