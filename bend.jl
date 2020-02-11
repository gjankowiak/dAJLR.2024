module Bend

import LinearAlgebra
const LA = LinearAlgebra

import SparseArrays
const SA = SparseArrays

import JankoUtils: spdiagm_const

const SubA = SubArray{Float64,1,Array{Float64,1}}

struct Params
    N::Int64
    M::Float64
    epsilon::Float64
    beta::Function
    beta_prime::Function
    beta_second::Function
end

mutable struct Candidate
    ρ::SubA
    θ::SubA
    λ1::Float64
    λ2::Float64
    λ3::Float64
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

function X2candidate(P::Params, X)
    return Candidate(view(X, 1:P.N),
                     view(X, P.N+1:2*P.N-1)
                     X[2*P.N],
                     X[2*P.N+1],
                     X[2*P.N+2])
end

function candidate2X(c::Candidate)
    return c.parent
end

"""
    Structure of X:
    * ρ = X[1:P.N]
    * θ = X[P.N+1, 2*P.N-1]
    * λ1 = X[2*P.N]
    * λ2 = X[2*P.N+1]
    * λ3 = X[2*P.N+2]
"""

function minimize_energy(P::Params, Xinit::Vector{Float64}; tol::Float64=1e-8, relaxation::Float64=1.0, max_iter::Int64=100)
    # One time initialization of finite differences matrices
    matrices = assemble_fd_matrices(P)

    # Initialization
    err = Inf
    X = copy(Xinit)
    n = 0

    b = compute_residuals(P, matrices, X)

    # Loop until tolerance or max. iterations is met
    while (err > tol) && (n < max_iter)

        # Assemble
        A = assemble_inner_system(P, matrices, X)

        # Solve
        δX = A\b

        # Update
        X += relaxation*δX
        b = compute_residuals(P, matrices, X)

        # Compute residual norm
        err = LA.norm(b)
        n += 1
    end
    return (X, err)
end

function assemble_inner_system(P::Params, matrices::FDMatrices, X::Vector{Float64}, relaxation)
    N = P.N

    # Easier access to ρ, θ, λi
    c = X2candidate(P, X)

    # Compute beta and derivatves
    beta_prime_ρ = P.beta_prime(c.ρ)
    beta_second_ρ = P.beta_second(c.ρ)
    beta_phalf = P.beta(matrices.M_phalf*c.ρ)[2:end]
    beta_mhalf = P.beta(matrices.M_mhalf*c.ρ)[2:end]
    beta_prime_phalf = P.beta_prime(matrices.M_phalf*c.ρ)
    beta_prime_mhalf = P.beta_prime(matrices.M_mhalf*c.ρ)

    # Compute upstream and downstream finite differences of θ
    θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, θ)
    θ_prime_centered = compute_centered_fd_θ(P, matrices, c.θ)

    A_ρρ = P.epsilon^2 * matrices.D2 - 0.5*beta_second_ρ.*θ_prime_centered.^2*Matrix{Float64}(LA.I, N, N)
    A_ρθ = 2*beta_prime_ρ.*θ_prime_centered*matrices.D1c
    A_ρλ3 = -ones(N)

    A_θρ  = ((beta_prime_phalf.*θ_prime_down.*matrices.M_phalf - beta_prime_mhalf.*θ_prime_up.*matrices.M_mhalf)[2:N,1:N])/Δs
    A_θθ  = (beta_phalf.*matrices.D1[2:N,:] - beta_mhalf.*matrices.D1[1:N-1,:])/Δs + (c.λ1*cos.(c.θ) - c.λ2*sin.(c.θ))*Matrix{Float64}(LA.I, N-1, N-1)
    A_θλ1 = sin.(c.θ)*ones(N-1)
    A_θλ2 = -cos.(c.θ)*ones(N-1)

    A_c1θ = -P.Δs*sin.(c.θ)
    A_c2θ = P.Δs*cos.(c.θ)
    A_c3ρ = P.Δs

    A = [[ A_ρρ     A_ρθ  zeros(N) zeros(N)      A_ρλ3 ];
         [ A_θρ     A_θθ      A_λ1     A_λ2 zeros(N-1) ];
         [ zeros(N) A_c1θ        0        0          0 ];
         [ zeros(N) A_c2θ        0        0          0 ];
         [ A_c3ρ    zeros(N-1)   0        0          0 ]]

    return A
end

function compute_residuals(P::Params, matrices::FDMatrices, X::Vector{Float64})
    # Easier access to ρ, θ, λi
    c = X2candidate(P, X)

    # Compute beta and derivatves
    beta_prime_ρ = P.beta_prime(c.ρ)
    beta_phalf = P.beta(matrices.M_phalf*c.ρ)[2:N]
    beta_mhalf = P.beta(matrices.M_mhalf*c.ρ)[2:N]

    # Compute upstream and downstream finite differences of θ
    θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, θ)

    # Compute residuals for ρ
    b_ρ = P.epsilon^2 * (matrices.D2 * c.ρ) - 0.5*beta_prime_ρ.*(compute_centered_fd_(P, matrices, c.θ)).^2 - P.λ3*ones(P.N)
    # Compute residuals for θ
    b_θ = (beta_phalf.*θ_prime_down - beta_mhalf.*θ_prime_up)/P.Δs + c.λ1*sin.(c.θ) - c.λ2*cos(θ)

    # Assemble with constraint residuals
    return [b_ρ; b_θ, P.Δs*sum(c.ρ) - P.M, P.Δs*sum(cos.(c.θ)), P.Δs*sum(sin.(c.θ))]
end

function compute_fd_θ(P::Params, matrices::FDMatrices, θ::SubA)
    fd = matrices.D1*θ + matrics.D1_rhs
    return (view(fd, 1:N-1), view(fd, 2:N))
end

function compute_centered_fd_θ(P::Params, matrices::FDMatrices, θ::SubA)
    return matrices.D1c*θ + matrics.D1c_rhs
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
                      [zeros(N-1); 1/Δs],

                      # Matrix for centered finite differences
                      (SA.spdiagm(-2  => -ones(N-2),
                                  -1  => zeros(N-1),
                                  0   => ones(N-1),
                                  N-2 => [-1])[1:N,1:N-1])*0.5/Δs,

                      # Affine term for centered finite differences
                      [π/Δs; zero(N-2); π/Δs]

                      # Matrix for 2nd order finite difference
                      spdiagm_const([1.0, -2.0, 1.0], [-1, 0, 1], N)/Δs^2
                     )
end

end # module
