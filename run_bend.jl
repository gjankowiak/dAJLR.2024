push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot

function beta_1(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return @. rho^2+1
end

function beta_1_prime(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return 2*rho
end

function beta_1_second(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return 2*ones(size(rho))
end

function beta_2(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return @. 3-rho^2
end

function beta_2_prime(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return -2*rho
end

function beta_2_second(rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return -2*ones(size(rho))
end

N = 250
Δs = 2π/N
M = 2π
epsilon = 1.0
epsilon = 5e-2

tol = 1e-7
max_iter = 8000
rel = 1e-3

# P_1 = Bend.Params(N, Δs, M, epsilon, beta_1, beta_1_prime, beta_1_second)
# P_2 = Bend.Params(N, Δs, M, epsilon, beta_2, beta_2_prime, beta_2_second)

beta0 = 1
rho0  = M/2π
m     = 1
h     = 1
j     = 4
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

print("Critical epsilon:")
println(eps_c)

# beta, beta_prime, beta_second = quadratic_beta(beta0, rho0, m, h)
P_3 = Bend.Params(N, Δs, M, epsilon, beta, beta_prime, beta_second)

P = P_3

Plotting.init()

# Circle
# Xcircle = Bend.initial_data(P, 1, 1)
# Plotting.plot(P, Xcircle, label="circle")

# Ellipsis
Xinit = Bend.initial_data(P, 1, 1, 3, false)

# Plotting.plot(P, Xinit, label="init")
# Plotting.s()

# display(Xinit)
println(P)

res = Bend.minimize_energy(P, Xinit; tol=tol, max_iter=max_iter, relaxation=rel)

print("Tolerance met: ")
println(res.converged)

Plotting.plot_result(P, res, label="solution")
Plotting.s()

# minimizor = Bend.minimizor(P, Xinit; max_iter=10, relaxation=1e-1)

# X, err = Xinit, Inf

# while true
    # (X, b, err, done) = minimizor()
    # if done
        # break
    # end
    # c = Bend.X2candidate(P, X)
    # println("θ")
    # display(c.θ)
    # println()
    # println(b[end-2])
    # Plotting.plot(P, X)
    # Plotting.s()
# end

# display(X)
# println()
# println(err)
