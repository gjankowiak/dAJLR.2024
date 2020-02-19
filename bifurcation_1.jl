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
Xinit = Bend.initial_data(P, 1, 1, 5)
Xcur = copy(Xinit)

print("Critical epsilon:")
println(eps_c)

delta_eps = 10. .^range(-5, -1, length=10)
epss = eps_c .+ [0.1; 0; -delta_eps]
neps = length(epss)

idx_first_bif = findfirst(x -> x == eps_c, epss)

amp_ρ = zeros(neps)
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
    (X, err) = Bend.minimize_energy(P, Xcur; tol=tol, max_iter=max_iter, relaxation=rel)
    c = Bend.X2candidate(P, X)
    extrema_rho = extrema(c.ρ)
    amp_ρ[i] = extrema_rho[2] - extrema_rho[1]
    errs[i] = err
    Xcur .= X
    println(err)
    # Plotting.plot(P, X)
end

data = [epss amp_ρ errs]

DelimitedFiles.writedlm("bifurcation_1_amplitude.txt", data, ",")
