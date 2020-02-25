push!(LOAD_PATH, "src")

import Bend
import DelimitedFiles
# import Plotting

N = 250
Δs = 2π/N
M = 2π
epsilon = 1.0
epsilon = 0.2

tol = 1e-7
max_iter = 100000
rel = 2e-2

beta0 = 1
rho0  = M/2π
m     = 1
h     = 1
j     = 2
eps_c = sqrt((m^2 - h/2)/j^2)

function beta(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return beta0*(1 .+ m*(rho.-rho0) + 0.5*h*(rho.-rho0).^2)
end

function beta_prime(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return beta0*(m*ones(size(rho)) + h*(rho.-rho0))
end

function beta_second(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return beta0*h*ones(size(rho))
end

P_3 = Bend.Params(N, Δs, M, epsilon, beta, beta_prime, beta_second)

P = P_3

# Ellipsis
# Xinit = Bend.initial_data(P, 1, 1, 5)
Xinit = Bend.initial_data(P, 1.1, 1, 1, false)
Xcur = copy(Xinit)

print("Critical epsilon:")
println(eps_c)

delta_eps = [10. .^range(-5, -1, length=10); collect(0.1.+0.05*(1:4))]
epss = eps_c .+ [0.1; 0; -delta_eps; ]
neps = length(epss)

idx_first_bif = findfirst(x -> x == eps_c, epss)

min_ρ = zeros(neps)
max_ρ = zeros(neps)
energies = zeros(neps)
errs = zeros(neps)

for i in 1:neps
    global X
    P.epsilon = epss[i]
    print(i)
    print(" ")
    println(P.epsilon)
    if i == idx_first_bif+1
        Xcur .= Xinit
    end
    res = Bend.minimize_energy(P, Xcur; atol=1e-10, max_iter=max_iter, relaxation=rel, adapt=true)
    X = res.sol
    iter = res.iter
    c = Bend.X2candidate(P, X)
    min_ρ[i], max_ρ[i] = extrema(c.ρ)
    energies[i] = res.energy_i[iter]
    errs[i] = res.residual_i[iter]
    Xcur .= X
    print("Residual: "); println(errs[i])
    print("Iterations: "); println(res.iter)
    print("Tolerance reached: "); println(res.converged)
    # Plotting.plot(P, X)
end

data = [epss min_ρ max_ρ energies errs]

DelimitedFiles.writedlm("bifurcation_1_amplitude.txt", data, ",")
