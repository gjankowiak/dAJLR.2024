push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot
import Serialization
import Arpack


N = 2*2*3*5*7
N = 2*2*2*2*2*2*7
print("N: ")
println(N)
Δs = 2π/N
M = 2π

epsilon = 2e-2
potential_range = 0.0
rho_max = -1e-12

center_rho = false

atol = 1e-8
rtol = 1e-13
max_iter = 3000
step_size = 1e-1

beta_0 = 1
beta_rho0  = M/2π
beta_m     = 0.4
beta_h     = 0.3
beta_k     = 0
beta_j     = 2

bigZ(m,h) = -m^6 + 8*m^4*h - 4*m^2*h^2 + h^3/2

bigZ2(m, h) = h^3/2 - 9*h^2*m^2 + 18*h*m^4 - 7*m^6

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
eps_delta = 1e-3
Xinit = Bend.initial_data(P, 1, 1, pulse=2, pulse_amplitude=7*sqrt(eps_delta), reverse_phase=false, only_rho=false)
# Xinit = Bend.initial_data_smooth(P, sides=8, smoothing=0.4, reverse_phase=true, only_rho=true)
#
# Xinit = Serialization.deserialize("case_6_branch_1_new.dat.bak")[1]
# Xinit = Serialization.deserialize("case_8_branch_1.dat")[end]
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
epsilons = [0.136]
# epsilons = collect(range(0.136, 0.135, step=-0.0001))
# epsilons = collect(range(0.135, 0.13, step=-0.001))
# epsilons = collect(range(0.13, 0.05, step=-0.01))
epsilons = collect(range(0.05, 0.01, step=-0.005))

# case 7
# # branch 1
epsilons = collect(range(0.62, 0.72, step=0.01))
epsilons = collect(range(0.73, 0.82, step=0.01))
epsilons = collect(range(0.87, 1.2, step=0.05))
# # branch 2
epsilons = [0.40825]
epsilons = collect(range(0.4083, 0.409, step=0.0001))
epsilons = collect(range(0.41, 0.45, step=0.0025))
epsilons = collect(range(0.46, 0.66, step=0.01))

# case 8
epsilon0 = 0.612372435695
epsilons = [epsilon0 + eps_delta]

# case m = 0.4, h = 0.3 # subcritical
P9 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 0.4, 0.3, 0, 2,
                 potential_range, center_rho)


# case m = 0.4, h = 0.04 # subcritical
P10 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 0.4, 0.04, 0, 2,
                 potential_range, center_rho)

# case m = 0.4, h = 0.2 # supercritical
P11 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 0.4, 0.2, 0, 2,
                 potential_range, center_rho)

P12 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1, 1.7, 0, 2,
                 potential_range, center_rho)

P13 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1, 1.71, 0, 2,
                 potential_range, center_rho)

P14 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1, 0.515, 0, 2,
                 potential_range, center_rho)
# P14 = Bend.Params(N, Δs,
                 # M, epsilon, rho_max,
                 # beta_0, beta_rho0, 1, 0.505, 0, 2,
                 # potential_range, center_rho)

P15 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1, 0.525, 0, 2,
                 potential_range, center_rho)

P16 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1.9, 3, 0, 2,
                 potential_range, center_rho)

P17 = Bend.Params(N, Δs,
                 M, epsilon, rho_max,
                 beta_0, beta_rho0, 1.9, 3, 0, 3,
                 potential_range, center_rho)

P16_bounded = Bend.Params(N, Δs,
                 M, epsilon, 4.1,
                 beta_0, beta_rho0, 1.9, 3, 0, 2,
                 1e-1, center_rho)

P17_bounded = Bend.Params(N, Δs,
                 M, epsilon, 4.1,
                 beta_0, beta_rho0, 1.9, 3, 0, 3,
                 1e-1, center_rho)


function iterate_filename(pattern)
    for i in 1:1000
        s, s_P = string(replace(pattern, "#" => i), ".dat"), string(replace(pattern, "#" => i), "_P.dat")
        if (isfile(s) || isfile(s_P))
            continue
        end
        return s, s_P
    end
end

function follow_branch(P, eps_delta, n_samples; eps_start=0.0, eps_direction=0.0, prefix="result_#", fn_xinit="")
    Xs = []
    Ps = []

    Z = bigZ2(P.beta_m, P.beta_h)
    eps_c = sqrt((P.beta_m^2 - P.beta_h/2)/P.beta_j^2)
    amp = sqrt(abs(P.beta_j^2*4*(2*P.beta_m^2 - P.beta_h)/Z))

    if eps_start > 0.0
        if eps_direction == 0.0
            throw("Please provide eps_direction")
        end
        epsilons = collect(eps_start .+ sign(eps_direction).*eps_delta*(0:n_samples-1))
    else
        epsilons = collect(eps_c .- sign(Z).*eps_delta*(1:n_samples))
    end

    @show epsilons


    @show Z
    println("j:", P.beta_j, ", eps: ", eps_c)
    println(", amplitude: ", amp)

    if fn_xinit == ""
        Xinit = Bend.initial_data(P, 1, 1, pulse=P.beta_j, pulse_amplitude=1e0*amp*sqrt(abs(eps_c - epsilons[1])), reverse_phase=false, only_rho=false)
    else
        Xinit = Serialization.deserialize(fn_xinit)[end]
        @show Xinit
    end
    Xx = Base.copy(Xinit)

    Xx = Base.copy(Xinit)

    f = Plotting.init(P)
    Plotting.plot(f, P, Xx)

    for (k,e) in enumerate(epsilons)

        if e == epsilons[1]
            Xx .= Xinit
        end

        P.epsilon = e
        println(k, "/", n_samples, " - epsilon: ", P.epsilon)

        minimizor = Bend.minimizor(P, Xx, solver_params)

        local res

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

            if res.finished
                print("Finished, converged: ")
                println(res.converged)
                print("Energy: ")
                print(res.energy_i[n])
                print(", residual norm: ")
                print(res.residual_norm_i[n])
                print(", min/max ρ: ")
                println(extrema(view(Xx, 1:P.N)))
                push!(Xs, res.sol)
                push!(Ps, Bend.copy(P))
                break
            end
        end
    end

    s, s_P = iterate_filename(prefix)

    Serialization.serialize(s, Xs)
    Serialization.serialize(s_P, Ps)
end

function follow_branch_decoupled(P, eps_delta, n_samples; eps_start=0.0, eps_direction=0.0, prefix="result_#", fn_xinit="")
    Xs = []
    Ps = []

    if P.beta_m != 0 || P.beta_h >= 0
        throw("Parameter m is not zero or h is non negative")
    end
    eps_c = sqrt(- P.beta_h/2/P.beta_j^2)
    amp = sqrt(8*P.beta_j^2/P.beta_h^2)

    if eps_start > 0.0
        if eps_direction == 0.0
            throw("Please provide eps_direction")
        end
        epsilons = collect(eps_start .+ sign(eps_direction).*eps_delta*(0:n_samples-1))
    else
        epsilons = collect(eps_c .+ eps_delta*(1:n_samples))
    end

    @show epsilons


    println("j:", P.beta_j, ", eps: ", eps_c)
    println(", amplitude: ", amp)

    if fn_xinit == ""
        Xinit = Bend.initial_data(P, 1, 1, pulse=P.beta_j, pulse_amplitude=amp*sqrt(abs(eps_c - epsilons[1])), reverse_phase=false, only_rho=false)
    else
        Xinit = Serialization.deserialize(fn_xinit)[end]
    end
    Xx = Base.copy(Xinit)

    f = Plotting.init(P)
    Plotting.plot(f, P, Xx)

    for (k,e) in enumerate(epsilons)

        if e == epsilons[1]
            Xx .= Xinit
        end

        P.epsilon = e
        println(k, "/", n_samples, " - epsilon: ", P.epsilon)

        minimizor = Bend.minimizor(P, Xx, solver_params)

        local res

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

            if res.finished
                print("Finished, converged: ")
                println(res.converged)
                push!(Xs, res.sol)
                push!(Ps, Bend.copy(P))
                break
            end
        end
    end

    s, s_P = iterate_filename(prefix)

    Serialization.serialize(s, Xs)
    Serialization.serialize(s_P, Ps)

end

function find_critical_point(P, eps_delta)
    Xs = []
    Ps = []

    eps_c = sqrt((P.beta_m^2 - P.beta_h/2)/P.beta_j^2)
    Z = bigZ2(P.beta_m, P.beta_h)
    amp = sqrt(abs(P.beta_j^2*4*(2*P.beta_m^2 - P.beta_h)/Z))

    println("2m^2 - h = ", 2*P.beta_m^2 - P.beta_h)

    print("j:", P.beta_j, ", eps: ", eps_c)
    print(", amplitude: ", amp)
    if Z > 0
        println(", supercritical")
    else
        println(", subcritical")
    end

    Xinit = Bend.initial_data(P, 1, 1, pulse=2, pulse_amplitude=0.5amp*sqrt(abs(eps_delta)), reverse_phase=false, only_rho=false)
    Xx = Base.copy(Xinit)

    f = Plotting.init(P)
    Plotting.plot(f, P, Xx)

    P.epsilon = eps_c - sign(Z)*eps_delta
    P.epsilon = eps_c - eps_delta
    println(P.epsilon)

    minimizor = Bend.minimizor(P, Xx, solver_params)

    local res

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
            push!(Xs, res.sol)
            push!(Ps, Bend.copy(P))
            # if e == epsilons[1]
                # Serialization.serialize("xstart_2.dat", res.sol)
            # end
            break
        end
    end
end

# find_critical_point(P13, 1e-5)
eps_delta = 1e-2
# follow_branch(P16, 2.5e-2, 8, eps_start=0.22629-2.5e-2, eps_direction=-1)
eps_delta = 1e-2
# follow_branch(P17, eps_delta, 8, eps_start=0.109194-eps_delta, eps_direction=-1)

eps_delta = 5e-4
# follow_branch(P16_bounded, eps_delta, 10)

eps_start = 0.7212919523166975
eps_delta = 5e-3

eps_start = 0.6712919523166975
eps_delta = 6e-2

eps_start = 0.07129195231669753
eps_delta = 5e-3

# follow_branch(P16_bounded, eps_delta, 10, eps_start=eps_start-eps_delta, eps_direction=-1)

P_decoupled = Bend.Params(N, Δs,
                 M, epsilon, 4.1,
                 beta_0, beta_rho0, 0, -1, 0, 1,
                 0.0, center_rho)
# follow_branch_decoupled(P_decoupled, 1e-4, 10, prefix="decoupled_j_1_#")

eps_start = 0.7081067811865476
eps_delta = 1e-3

# follow_branch_decoupled(P_decoupled, eps_delta, 10, eps_start=eps_start+eps_delta, eps_direction=1.0, prefix="decoupled_j_1_#")

eps_start = 0.7181067811865476
eps_delta = 1e-1

# eps_start = 1.6781067811
# eps_delta = 1e-2

eps_start = 3.2181067811865476
eps_delta = 5e-1

# follow_branch_decoupled(P_decoupled, eps_delta, 50, eps_start=eps_start, eps_direction=1.0, prefix="decoupled_j_1_#", fn_xinit="decoupled_j_1_5.dat")

P_decoupled = Bend.Params(N, Δs,
                 M, epsilon, 4.1,
                 beta_0, beta_rho0, 0, -1, 0, 2,
                 0.0, center_rho)

# Decoupled
# follow_branch_decoupled(P_decoupled, 1e-4, 10)

eps_start = 0.3545533905932738
eps_delta = 1e-3
# follow_branch_decoupled(P_decoupled, eps_delta, 10, eps_start=eps_start+eps_delta, eps_direction=1.0)

eps_start = 0.3645533905932738
eps_delta = 1e-2
# follow_branch_decoupled(P_decoupled, eps_delta, 10, eps_start=eps_start+eps_delta, eps_direction=1.0)

P_decoupled = Bend.Params(N, Δs,
                 M, epsilon, 4.1,
                 beta_0, beta_rho0, 0, -1, 0, 3,
                 0.0, center_rho)

# Decoupled
# follow_branch_decoupled(P_decoupled, 1e-4, 10)

eps_start = 0.23670226039551584
eps_delta = 1e-3
# follow_branch_decoupled(P_decoupled, eps_delta, 10, eps_start=eps_start+eps_delta, eps_direction=1.0)

# eps_start = 0.24670226039551585
# eps_delta = 1e-2
# follow_branch_decoupled(P_decoupled, eps_delta, 10, eps_start=eps_start+eps_delta, eps_direction=1.0)

eps_delta = 1e-2

Pm1m1 = Bend.Params(N, Δs,
                    M, epsilon, 0,
                    beta_0, beta_rho0, -1, -1, 0, 2,
                    0.0, center_rho)


eps_start = 0.5
eps_delta = 1e-4

follow_branch(Pm1m1, eps_delta, 5, prefix="Pm1m1_j_2_#")
# follow_branch(Pm1m1, eps_delta, eps_start=eps_start - eps_delta, eps_direction=-1, 5, prefix="Pm1m1_j_2_#")
# follow_branch(Pm1m1, eps_delta, eps_start=eps_start + eps_delta, eps_direction=1, 5, prefix="Pm1m1_j_2_#")


Pm2m2 = Bend.Params(N, Δs,
                    M, epsilon, 0,
                    beta_0, beta_rho0, -2, -2, 0, 2,
                    0.0, center_rho)

eps_delta = 1e-4
# follow_branch(Pm2m2, eps_delta, 25, prefix="Pm2m2_j_2_#")

eps_start = 1.1280339
eps_delta = 1e-3
#
# follow_branch(Pm2m2, eps_start=eps_start, eps_direction=-1, eps_delta, 10, prefix="Pm2m2_j_2_#")

eps_start = 1.368033988749895
eps_delta = 1e-1

# follow_branch(Pm2m2, eps_start=eps_start, eps_direction=1, eps_delta, 25, prefix="Pm2m2_j_2_#", fn_xinit="Pm2m2_j_2_1.dat")

Pm2m2 = Bend.Params(N, Δs,
                    M, epsilon, 0,
                    beta_0, beta_rho0, -2, -2, 0, 3,
                    0.0, center_rho)

# follow_branch(Pm1m1, eps_delta, eps_start=eps_start + eps_delta, eps_direction=1, 5, prefix="Pm1m1_j_2_#")

# follow_branch(Pm2m2, eps_delta, 25, prefix="Pm2m2_j_2_#")
eps_start = 0.7703559924999299
eps_delta = 1e-2
# follow_branch(Pm2m2, eps_start=eps_start, eps_direction=1, eps_delta, 25, prefix="Pm2m2_j_2_#", fn_xinit="Pm2m2_j_2_1.dat")

eps_start = 1.01035599249993
eps_delta = 1e-2
# follow_branch(Pm2m2, eps_start=eps_start, eps_direction=1, eps_delta, 25, prefix="Pm2m2_j_2_#", fn_xinit="Pm2m2_j_2_2.dat")

eps_start = 1.25035599249993
eps_delta = 3e-2
# follow_branch(Pm2m2, eps_start=eps_start, eps_direction=1, eps_delta, 25, prefix="Pm2m2_j_2_#", fn_xinit="Pm2m2_j_2_3.dat")

