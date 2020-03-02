module Bend

import LinearAlgebra
const LA = LinearAlgebra

import SparseArrays
const SA = SparseArrays

import JankoUtils: spdiagm_const
import EvenParam

const SubA = SubArray{Float64,1,Array{Float64,1}}

mutable struct Params
    # Discretization parameters
    N::Int64
    Δs::Float64

    # Model parameters
    M::Float64
    epsilon::Float64
    ρ_max::Float64

    # Beta parameters
    beta_0::Float64
    beta_rho0::Float64
    beta_m::Float64
    beta_h::Float64
    beta_k::Float64
    beta_j::Int64

    # rho bound parameters
    potential_range::Float64

    # solver parameters
    center_ρ::Bool
end

mutable struct SolverParams
    atol::Float64
    rtol::Float64
    max_iter::Int64
    step_size::Float64
    adapt::Bool
    min_step_size::Float64
    max_step_size::Float64
    step_down_threshold::Float64
    step_up_threshold::Float64
    step_factor::Float64
end

default_SP = SolverParams(1e-8,    # atol
                          1e-12,   # rtol
                          1000,    # max_iter
                          1e-3,    # step_size
                          false,   # adapt
                          1e-5,    # min_step_size
                          1e-1,    # max_step_size
                          1e-1,    # step_down_threshold
                          1e-4,    # step_up_threshold
                          1.2)     # step_factor

mutable struct Candidate
    ρ::Union{Vector{Float64},SubA}
    θ::Union{Vector{Float64},SubA}
    λx::Float64
    λy::Float64
    λM::Float64
    λcm::Float64
end

mutable struct History
    energy_prev::Float64
    residual_prev::Vector{Float64}
end

struct Result
    sol::Vector{Float64}
    iter::Int64
    energy_i::Vector{Float64}
    residual_norm_i::Vector{Float64}
    converged::Bool
    finished::Bool
end

mutable struct Relaxation
    r::Float64
    threshold::Float64
    min_r::Float64
    max_r::Float64
end

struct FDMatrices
    M_phalf::SA.SparseMatrixCSC{Float64}
    M_mhalf::SA.SparseMatrixCSC{Float64}
    D1::SA.SparseMatrixCSC{Float64}
    D1_rhs::Vector{Float64}
    D1c::SA.SparseMatrixCSC{Float64}
    D1c_rhs::Vector{Float64}
    D2::SA.SparseMatrixCSC{Float64}
end

function copy(P::Params)
    return Params(
    P.N, P.Δs, P.M,
    P.epsilon, P.ρ_max, P.beta_0,
    P.beta_rho0, P.beta_m, P.beta_h,
    P.beta_k, P.beta_j, P.potential_range,
    P.center_ρ)
end

function prompt_yes_no(s::String, default_yes::Bool=false)
    while true
        print(s)
        if default_yes
            print(" [Y/n]: ")
        else
            print(" [y/N]: ")
        end
        i = readline()
        if (default_yes && i == "") || (i == "y")
            return true
        elseif (!default_yes && i == "") || (i == "n")
            return false
        end
    end
end

function X2candidate(P::Params, X; copy::Bool=false)
    if P.center_ρ
        if copy
            return Candidate(X[1:P.N],
                             X[P.N+1:2*P.N-1],
                             X[2*P.N],
                             X[2*P.N+1],
                             X[2*P.N+2],
                             X[2*P.N+3])
        else
            return Candidate(view(X, 1:P.N),
                             view(X, P.N+1:2*P.N-1),
                             X[2*P.N],
                             X[2*P.N+1],
                             X[2*P.N+2],
                             X[2*P.N+3])
        end
    else
        if copy
            return Candidate(X[1:P.N],
                             X[P.N+1:2*P.N-1],
                             X[2*P.N],
                             X[2*P.N+1],
                             X[2*P.N+2],
                             0)
        else
            return Candidate(view(X, 1:P.N),
                             view(X, P.N+1:2*P.N-1),
                             X[2*P.N],
                             X[2*P.N+1],
                             X[2*P.N+2],
                             0)
        end
    end
end

function candidate2X(c::Candidate)
    return c.parent
end

function compute_beta(P::Params, rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return P.beta_0*(1 .+
                         P.beta_m*(rho.-P.beta_rho0) +
                     0.5*P.beta_h*(rho.-P.beta_rho0).^2 +
                     1/4*P.beta_k*(rho.-P.beta_rho0).^4)
end

function compute_beta_prime(P::Params, rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return P.beta_0*(P.beta_m*ones(size(rho)) +
                     P.beta_h*(rho.-P.beta_rho0) +
                     P.beta_k*(rho.-P.beta_rho0).^3)
end

function compute_beta_second(P::Params, rho::Union{Vector{Float64},SubArray{Float64,1,Array{Float64,1}}})
    return P.beta_0*(  P.beta_h*ones(size(rho)) +
                     3*P.beta_k*(rho.-P.beta_rho0).^2)
end

"""
    Structure of X:
    * ρ = X[1:P.N]
    * θ = X[P.N+1, 2*P.N-1]
    * λx = X[2*P.N]
    * λy = X[2*P.N+1]
    * λM = X[2*P.N+2]
"""

function adapt_step(P::Params, SP::SolverParams, matrices::FDMatrices,
                    X::Vector{Float64}, δX::Vector{Float64}, current_residual::Vector{Float64},
                    current_step_size::Float64, history::History)

    local energy, residual

    Xnew = zeros(size(X))
    t = current_step_size

    while true
        #
        # Update
        Xnew .= X + t*δX
        cnew = X2candidate(P, Xnew)

        # Check that rho is within admissible bounds
        # IF not, reduce step size
        if ((P.potential_range > 0) && (maximum(cnew.ρ) > P.ρ_max))
            t /= 2
            print("Rho above bound, reducing relaxation parameter to: ")
            println(t)
            if t < SP.min_step_size
                error("Relaxation parameter became too small")
            end
            continue
        end

        # Compute residual and energy
        residual = compute_residual(P, matrices, Xnew)
        energy = compute_energy(P, matrices, Xnew)

        if (energy - history.energy_prev > history.energy_prev*1e2)
            t /= SP.step_factor
            print("Reducing relaxation parameter to: ")
            println(t)
            if t < SP.min_step_size
                error("Relaxation parameter became too small")
            end
        else
            if ((0 < (history.energy_prev - energy)) &&
                (LA.norm(residual) > 1e-3) &&
                (LA.norm(residual - history.residual_prev)/LA.norm(residual) < SP.step_up_threshold) &&
                (t < SP.max_step_size))
                t = min(SP.step_factor*t, SP.max_step_size)
                print("Increasing relaxation parameter to: ")
                println(t)
            else
                # print("residual ratio: ")
                # println(LA.norm(residual - history.residual_prev)/LA.norm(residual))
            end
            @goto finish
        end
    end

    @label finish
    history.energy_prev = energy
    history.residual_prev = current_residual

    return Xnew, residual, t
end

function minimizor(P::Params, Xinit::Vector{Float64}, SP::SolverParams=default_SP)
    matrices = assemble_fd_matrices(P)

    # Initialization
    history = History(0, zeros(2*P.N+2 + P.center_ρ))

    X = Base.copy(Xinit)
    n = 1
    step_size = SP.step_size

    residual = compute_residual(P, matrices, X)

    energy_i = zeros(SP.max_iter)
    residual_norm_i = zeros(SP.max_iter)

    energy_i[1] = compute_energy(P, matrices, X)
    residual_norm_i[1] = LA.norm(residual)

    history.energy_prev = energy_i[1]
    history.residual_prev .= residual

    # Loop until tolerance or max. iterations is met
    function ()
        n += 1
        # Assemble
        A = assemble_inner_system(P, matrices, X)

        # Solve
        δX = A\-residual

        if SP.adapt
            Xnew, residual, step_size = adapt_step(P, SP, matrices,
                                                   X, δX,
                                                   residual, step_size,
                                                   history)
            X .= Xnew
        else
            # Update
            X += step_size*δX
            residual = compute_residual(P, matrices, X)
        end

        # Compute residual norm and energy
        residual_norm = LA.norm(residual)
        energy = compute_energy(P, matrices, X)
        energy_i[n] = energy
        residual_norm_i[n] = residual_norm

        converged = (LA.norm(residual - history.residual_prev) < SP.rtol) || (LA.norm(residual) < SP.atol)

        history.energy_prev = energy
        history.residual_prev = residual

        res = Result(X, n, energy_i, residual_norm_i, converged, converged || (n>=SP.max_iter))
        return res
    end
end

function compute_energy_split(P::Params, matrices::FDMatrices, X::Vector{Float64})
    c = X2candidate(P, X)
    ρ_dot = (circshift(c.ρ, -1) - c.ρ)/P.Δs
    θ_dot = compute_centered_fd_θ(P, matrices, c.θ)
    beta = compute_beta(P, c.ρ)
    return (0.5*P.epsilon^2*P.Δs*sum(ρ_dot.^2), 0.5P.Δs*sum(beta.*θ_dot.^2))
end

function compute_energy(P::Params, matrices::FDMatrices, X::Vector{Float64})
    c = X2candidate(P, X)
    ρ_dot = (circshift(c.ρ, -1) - c.ρ)/P.Δs
    θ_dot = compute_centered_fd_θ(P, matrices, c.θ)
    beta = compute_beta(P, c.ρ)
    E = 0.5*P.epsilon^2*P.Δs*sum(ρ_dot.^2) + 0.5P.Δs*sum(beta.*θ_dot.^2)
    if P.potential_range > 0
        (w, _, _) = compute_potential(P, X)
        E += P.Δs*sum(w)
    end
    return E
end

function assemble_inner_system(P::Params, matrices::FDMatrices, X::Vector{Float64})
    N = P.N
    Δs = P.Δs

    # Easier access to ρ, θ, λi
    c = X2candidate(P, X)

    # Compute beta and derivatves
    beta_prime_ρ = compute_beta_prime(P, c.ρ)
    beta_second_ρ = compute_beta_second(P, c.ρ)
    beta_phalf = compute_beta(P, matrices.M_phalf*c.ρ)[2:end]
    beta_mhalf = compute_beta(P, matrices.M_mhalf*c.ρ)[2:end]
    beta_prime_phalf = compute_beta_prime(P, matrices.M_phalf*c.ρ)
    beta_prime_mhalf = compute_beta_prime(P, matrices.M_mhalf*c.ρ)

    # Compute upstream and downstream finite differences of θ
    θ_prime_down_full, θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, c.θ)
    θ_prime_centered = compute_centered_fd_θ(P, matrices, c.θ)

    A_E1_ρ = P.epsilon^2 * matrices.D2 - 0.5*beta_second_ρ.*θ_prime_centered.^2 .*Matrix{Float64}(LA.I, N, N)
    if P.potential_range > 0
        (_, _, w_second) = compute_potential(P, X)
        A_E1_ρ .-= w_second .* Matrix{Float64}(LA.I, N, N)
    end
    # FIXME should be minus sign
    A_E1_θ = -beta_prime_ρ.*θ_prime_centered.*matrices.D1c

    A_E1_λM = ones(N)

    A_E2_ρ  = (beta_prime_phalf[2:N].*θ_prime_down.*matrices.M_phalf[2:N,:] - beta_prime_mhalf[2:N].*θ_prime_up.*matrices.M_mhalf[2:N,:])/Δs
    A_E2_θ  = (beta_phalf.*matrices.D1[2:N,:] - beta_mhalf.*matrices.D1[1:N-1,:])/Δs - (c.λx*cos.(c.θ) + c.λy*sin.(c.θ)).*Matrix{Float64}(LA.I, N-1, N-1)

    A_E2_λx = -sin.(c.θ)
    A_E2_λy =  cos.(c.θ)

    # FIXME
    # A_E2_ρ  *= 1e4
    # A_E2_θ  *= 1e4
    # A_E2_λx *= 1e4
    # A_E2_λy *= 1e4

    A_c1θ = -P.Δs*sin.(c.θ)'
    A_c2θ = P.Δs*cos.(c.θ)'
    A_c3ρ = P.Δs*ones(N)'

    if P.center_ρ
        A_E1_λcm = P.Δs.*(0:(N-1))
        A_c4ρ  = A_E1_λcm'
        A = [[ A_E1_ρ      A_E1_θ   zeros(N) zeros(N)      A_E1_λM     A_E1_λcm ];
             [ A_E2_ρ      A_E2_θ      A_E2_λx    A_E2_λy zeros(N-1) zeros(N-1) ];
             [ zeros(N)' A_c1θ         0        0          0          0 ];
             [ zeros(N)' A_c2θ         0        0          0          0 ];
             [ A_c3ρ     zeros(N-1)'   0        0          0          0 ];
             [ A_c4ρ     zeros(N-1)'   0        0          0          0 ]]
    else
        A = [[ A_E1_ρ    A_E1_θ        zeros(N) zeros(N)   A_E1_λM ];
             [ A_E2_ρ    A_E2_θ        A_E2_λx  A_E2_λy    zeros(N-1) ];
             [ zeros(N)' A_c1θ         0        0          0 ];
             [ zeros(N)' A_c2θ         0        0          0 ];
             [ A_c3ρ     zeros(N-1)'   0        0          0 ]]
    end

    return A
end

function compute_residual(P::Params, matrices::FDMatrices, X::Vector{Float64})
    N = P.N
    Δs = P.Δs

    # Easier access to ρ, θ, λi
    c = X2candidate(P, X)

    # Compute beta and derivatves
    beta_prime_ρ = compute_beta_prime(P, c.ρ)
    beta_phalf = compute_beta(P, matrices.M_phalf*c.ρ)[2:N]
    beta_mhalf = compute_beta(P, matrices.M_mhalf*c.ρ)[2:N]

    # Compute upstream and downstream finite differences of θ
    θ_prime_down_full, θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, c.θ)

    # Compute residuals for ρ
    b_E1 = P.epsilon^2 * (matrices.D2 * c.ρ) - 0.5*beta_prime_ρ.*(compute_centered_fd_θ(P, matrices, c.θ)).^2 + c.λM*ones(N)
    if P.potential_range > 0
        (_, w_prime, _) = compute_potential(P, X)
        b_E1 -= w_prime
    end
    # Compute residuals for θ
    b_E2 = (beta_phalf.*θ_prime_down - beta_mhalf.*θ_prime_up)/Δs - c.λx*sin.(c.θ) + c.λy*cos.(c.θ)

    # Assemble with constraint residuals
    if P.center_ρ
        return [b_E1; b_E2; Δs*(1+sum(cos.(c.θ))); Δs*sum(sin.(c.θ)); (Δs*sum(c.ρ) - P.M); Δs*sum(c.ρ.*(0:(N-1))) - P.M/2π]
    else
        return [b_E1; b_E2; Δs*(1+sum(cos.(c.θ))); Δs*sum(sin.(c.θ)); Δs*sum(c.ρ) - P.M]
    end
end

function compute_fd_θ(P::Params, matrices::FDMatrices, θ::Union{Vector{Float64}, SubA})
    fd = matrices.D1*θ + matrices.D1_rhs
    return (fd, view(fd, 1:P.N-1), view(fd, 2:P.N))
end

function compute_centered_fd_θ(P::Params, matrices::FDMatrices, θ::Union{Vector{Float64}, SubA})
    return matrices.D1c*θ + matrices.D1c_rhs
end

function assemble_fd_matrices(P::Params)
    N = P.N
    Δs = P.Δs

    return FDMatrices(
                      # Matrices for i+1/2 and i-1/2 values
                      spdiagm_const([0.5, 0.5], [0, 1], N),
                      spdiagm_const([0.5, 0.5], [-1, 0], N),

                      # Matrix for upstream and downstream 1st order finite diff
                      # Taking into account 0 and 2π boundary conditions
                      # Size (N x N-1)
                      # The finite differences can be computed using the
                      # compute_fd_θ function
                      (SA.spdiagm(-1 => -ones(N-1), 0 => ones(N-1))[1:N,1:N-1])/Δs,

                      # Affine term for up/downstream finite differences
                      [zeros(N-1); 2π/Δs],

                      # Matrix for centered finite differences
                      (SA.spdiagm(-2  => -ones(N-2),
                                  -1  => zeros(N-1),
                                  0   => ones(N-1),
                                  N-2 => [-1])[1:N,1:N-1])*0.5/Δs,

                      # Affine term for centered finite differences
                      [π/Δs; zeros(N-2); π/Δs],

                      # Matrix for 2nd order finite difference
                      spdiagm_const([1.0, -2.0, 1.0], [-1, 0, 1], N)/Δs^2
                     )
end

function initial_data_smooth(P::Params; sides::Int=1, smoothing::Float64, reverse_phase::Bool=false)
    if P.N % sides != 0
        error("N must be dividible by the number of sides")
    end

    if !(0 <= smoothing <= 1)
        error("smoothing parameter must be between 0 and 1")
    end

    k = Int64(P.N / sides)
    N_straight = floor(Int64, P.N/sides*(1-smoothing))

    thetas_1p = [zeros(N_straight); [2π/sides - i*P.Δs/smoothing for i in (k-N_straight-1):-1:0]]
    thetas = repeat(thetas_1p, sides) + repeat(2π/sides*(0:sides-1), inner=k)

    if reverse_phase
        rho = sin.(range(0, π, length=N_straight))
        rho /= (P.M/sides/(N_straight)) / sum(rho)
        rhos = repeat([rho; zeros(k-N_straight)], sides)
    else
        rhos = repeat([zeros(N_straight); P.M/(2π*smoothing)*ones(k-N_straight)], sides)
    end

    if P.center_ρ
        return [P.M/2π*ones(P.N); thetas[2:end]; 0; 0; 0; 0]
        # return [rhos; thetas[2:end]; 0; 0; 0; 0]
    else
        return [P.M/2π*ones(P.N); thetas[2:end]; 0; 0; 0]
        # return [rhos; thetas[2:end]; 0; 0; 0]
    end
end

function initial_data(P::Params, a::Real, b::Real; pulse::Int=1, pulse_amplitude::Real=2e-2, poly::Bool=false, reverse_phase::Bool=false)
    N = P.N

    # Construct the ellipsis and reparameterize it
    if poly && (pulse > 2)
        x = [real(exp(2π*im*k/pulse)) for k in 0:(pulse-1)]
        y = [imag(exp(2π*im*k/pulse)) for k in 0:(pulse-1)]
        xy = [x y]
    else
        t = collect(range(-π/2, 3π/2, length=10N+1)[1:10N])
        xy = [a*cos.(t) b*sin.(t)]
    end
    xy_even = EvenParam.reparam(xy; closed=true, new_N = N)

    # Recenter
    xy_even = xy_even .- [xy_even[1,1] xy_even[1,2]]

    # Rotate so that the first segment is parallel to the x axis
    xy_even_c = xy_even[:,1] + xy_even[:,2]*im
    alpha = angle(xy_even_c[2])
    xy_even_c = xy_even_c.*exp.(-im*alpha)

    # Compute tangential angles:w
    xy_even = [real.(xy_even_c) imag.(xy_even_c)]
    xy_diff = circshift(xy_even, -1) - xy_even
    xy_angles = angle.(xy_diff[:,1]+xy_diff[:,2]*im)

    thetas = zeros(N-1)
    thetas[1] = xy_angles[2]
    for i in 2:N-1
        thetas[i] = thetas[i-1] + rem(xy_angles[i+1] - thetas[i-1], Float64(π), RoundNearest)
    end

    t = collect(range(0, 2π, length=N+1)[1:N])
    if pulse > 0
        thetas += pulse_amplitude/pulse*sin.(pulse*t[2:N])
    end
    rhos = P.M/2π*ones(N)
    if pulse > 0
        if reverse_phase
            f = sin
        else
            f = cos
        end
        rhos -= pulse_amplitude*P.M/2π/pulse*f.(pulse*t)
    end

    if P.center_ρ
        return [rhos; thetas; 0; 0; 0; 0]
    else
        return [rhos; thetas; 0; 0; 0]
    end
end

function g(x::Vector{Float64}, α::Float64)
    return @. -min(α*x-1, 0.0)^2 * log(α*x)
end

function g_p(x::Vector{Float64}, α::Float64)
    return @. -(2α*min(α*x-1, 0.0)*log(α*x) + min(α*x-1, 0.0)^2/x)
end

function g_pp(x::Vector{Float64}, α::Float64)
    return @. -(2α^2*(x<(1/α))*log(α*x) + 4α*min(α*x-1, 0.0)/x - min(α*x-1, 0.0)^2/x^2)
end

function compute_potential(P::Params, X::Vector{Float64})
    α = P.potential_range
    ρ_max = P.ρ_max
    c = X2candidate(P, X)

    factor = 1e1

    @show minimum(ρ_max .- c.ρ)

    w = factor*g(ρ_max .- c.ρ, α)
    w_prime = factor*g_p(ρ_max .- c.ρ, α)
    w_second = factor*g_pp(ρ_max .- c.ρ, α)

    return (w, w_prime, w_second)
end

end # module
