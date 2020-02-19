module Bend

import LinearAlgebra
const LA = LinearAlgebra

import SparseArrays
const SA = SparseArrays

import JankoUtils: spdiagm_const
import EvenParam

const SubA = SubArray{Float64,1,Array{Float64,1}}

mutable struct Params
    N::Int64
    Δs::Float64
    M::Float64
    epsilon::Float64
    beta::Function
    beta_prime::Function
    beta_second::Function
end

mutable struct Candidate
    ρ::Union{Vector{Float64},SubA}
    θ::Union{Vector{Float64},SubA}
    λx::Float64
    λy::Float64
    λM::Float64
end

struct Result
    sol::Vector{Float64}
    iter::Int64
    energy_i::Vector{Float64}
    residual_i::Vector{Float64}
    converged::Bool
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

function X2candidate(P::Params, X; copy::Bool=false)
    if copy
        return Candidate(X[1:P.N],
                         X[P.N+1:2*P.N-1],
                         X[2*P.N],
                         X[2*P.N+1],
                         X[2*P.N+2])
    else
        return Candidate(view(X, 1:P.N),
                         view(X, P.N+1:2*P.N-1),
                         X[2*P.N],
                         X[2*P.N+1],
                         X[2*P.N+2])
    end
end

function candidate2X(c::Candidate)
    return c.parent
end

"""
    Structure of X:
    * ρ = X[1:P.N]
    * θ = X[P.N+1, 2*P.N-1]
    * λx = X[2*P.N]
    * λy = X[2*P.N+1]
    * λM = X[2*P.N+2]
"""

function minimizor(P::Params, Xinit::Vector{Float64}; tol::Float64=1e-8, relaxation::Float64=1.0, max_iter::Int64=100, dirichlet_rho::Bool=false)
    # One time initialization of finite differences matrices
    matrices = assemble_fd_matrices(P)

    # Initialization
    err = Inf
    X = copy(Xinit)
    n = 0

    b = compute_residuals(P, matrices, X)

    # Loop until tolerance or max. iterations is met
    function()
        # Assemble
        A = assemble_inner_system(P, matrices, X)
        if dirichlet_rho
            A[1,1] = 1e300
            b[1] = 0
        end

        # Solve
        δX = A\(-b)

        # Update
        X += relaxation*δX
        b = compute_residuals(P, matrices, X)

        # Compute residual norm
        err = LA.norm(b)

        n += 1

        if (err > tol) && (n < max_iter)
            (X, b, err, false)
        else
            (X, b, err, true)
        end
    end
end

function minimize_energy(P::Params, Xinit::Vector{Float64}; tol::Float64=1e-8, relaxation::Float64=1.0, max_iter::Int64=100, dirichlet_rho::Bool=false)
    # One time initialization of finite differences matrices
    matrices = assemble_fd_matrices(P)

    # Initialization
    err = Inf
    X = copy(Xinit)
    n = 0

    rel = relaxation
    b = compute_residuals(P, matrices, X)

    energy_i = zeros(max_iter)
    residual_i = zeros(max_iter)

    energy_i[1] = energy(P, matrices, X)
    residual_i[1] = LA.norm(b)

    # Loop until tolerance or max. iterations is met
    while (err > tol) && (n < max_iter)
        # Assemble
        A = assemble_inner_system(P, matrices, X)
        if dirichlet_rho
            A[1,1] = 1e300
            b[1] = 0
        end

        # Solve
        δX = A\-b

        if n == 2000
            rel *= 4
        end

        # Update
        X += rel*δX
        b = compute_residuals(P, matrices, X)

        # Compute residual norm
        err = LA.norm(b)

        n += 1

        energy_i[n] = energy(P, matrices, X)
        residual_i[n] = err
    end
    res = Result(X, n, energy_i, residual_i, n<max_iter)
    return res
end

function energy(P::Params, matrices::FDMatrices, X::Vector{Float64})
    c = X2candidate(P, X)
    ρ_dot = (circshift(c.ρ, -1) - c.ρ)/P.Δs
    θ_dot = compute_centered_fd_θ(P, matrices, c.θ)
    beta = P.beta(c.ρ)
    return 0.5*P.epsilon^2*P.Δs*sum(ρ_dot.^2) + 0.5P.Δs*sum(beta.*θ_dot.^2)
end

function assemble_inner_system(P::Params, matrices::FDMatrices, X::Vector{Float64})
    N = P.N
    Δs = P.Δs

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
    θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, c.θ)
    θ_prime_centered = compute_centered_fd_θ(P, matrices, c.θ)

    A_ρρ = P.epsilon^2 * matrices.D2 - 0.5*beta_second_ρ.*θ_prime_centered.^2 .*Matrix{Float64}(LA.I, N, N)
    A_ρθ = 2*beta_prime_ρ.*θ_prime_centered.*matrices.D1c
    # Δs factor for symmetry
    A_ρλM = -P.Δs*ones(N)

    A_θρ  = (beta_prime_phalf[2:N].*θ_prime_down.*matrices.M_phalf[2:N,:] - beta_prime_mhalf[1:N-1].*θ_prime_up.*matrices.M_mhalf[1:N-1,:])/Δs
    A_θθ  = (beta_phalf.*matrices.D1[2:N,:] - beta_mhalf.*matrices.D1[1:N-1,:])/Δs + (c.λx*cos.(c.θ) - c.λy*sin.(c.θ)).*Matrix{Float64}(LA.I, N-1, N-1)

    A_θλx = P.Δs*sin.(c.θ).*ones(N-1)
    A_θλy = -P.Δs*cos.(c.θ).*ones(N-1)

    # Δs factor for symmetry
    A_c1θ = -P.Δs*sin.(c.θ)'
    A_c2θ = P.Δs*cos.(c.θ)'
    A_c3ρ = P.Δs*ones(N)'

    A = [[ A_ρρ      A_ρθ   zeros(N) zeros(N)      A_ρλM ];
         [ A_θρ      A_θθ      A_θλx    A_θλy zeros(N-1) ];
         [ zeros(N)' A_c1θ         0        0          0 ];
         [ zeros(N)' A_c2θ         0        0          0 ];
         [ A_c3ρ     zeros(N-1)'   0        0          0 ]]

    # FIXME
    # M_nn = Matrix{Float64}(LA.I, N, N)

    # A = [[ M_nn         SA.spzeros(N,N-1)   SA.spzeros(N,1) SA.spzeros(N,1)  SA.spzeros(N,1) ];
         # [ SA.spzeros(N-1,N)      A_θθ      A_θλx    A_θλy SA.spzeros(N-1,1) ];
         # [ SA.spzeros(1,N)    A_c1θ         0        0          0 ];
         # [ SA.spzeros(1,N)    A_c2θ         0        0          0 ];
         # [ SA.spzeros(1,N)    SA.spzeros(1,N-1)   0        0          1 ]]

    return A
end

function compute_residuals(P::Params, matrices::FDMatrices, X::Vector{Float64})
    N = P.N
    Δs = P.Δs

    # Easier access to ρ, θ, λi
    c = X2candidate(P, X)

    # Compute beta and derivatves
    beta_prime_ρ = P.beta_prime(c.ρ)
    beta_phalf = P.beta(matrices.M_phalf*c.ρ)[2:N]
    beta_mhalf = P.beta(matrices.M_mhalf*c.ρ)[2:N]

    # Compute upstream and downstream finite differences of θ
    θ_prime_up, θ_prime_down = compute_fd_θ(P, matrices, c.θ)

    # Compute residuals for ρ
    b_ρ = P.epsilon^2 * (matrices.D2 * c.ρ) - 0.5*beta_prime_ρ.*(compute_centered_fd_θ(P, matrices, c.θ)).^2 - c.λM*ones(N)
    # Compute residuals for θ
    b_θ = (beta_phalf.*θ_prime_down - beta_mhalf.*θ_prime_up)/Δs + c.λx*sin.(c.θ) - c.λy*cos.(c.θ)

    # Assemble with constraint residuals
    return [b_ρ; b_θ; Δs*(1+sum(cos.(c.θ))); Δs*sum(sin.(c.θ)); Δs*sum(c.ρ) - P.M]
end

function compute_fd_θ(P::Params, matrices::FDMatrices, θ::SubA)
    fd = matrices.D1*θ + matrices.D1_rhs
    return (view(fd, 1:P.N-1), view(fd, 2:P.N))
end

function compute_centered_fd_θ(P::Params, matrices::FDMatrices, θ::SubA)
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

function initial_data(P::Params, a::Real, b::Real, pulse::Int=1, poly::Bool=false)
    N = P.N

    # Construct the ellipsis and reparameterize it
    t = collect(range(-π/2, 3π/2, length=N+1)[1:N])
    if pulse > 1
        if poly
            x = [real(exp(2π*im*k/pulse)) for k in 0:(pulse-1)]
            y = [imag(exp(2π*im*k/pulse)) for k in 0:(pulse-1)]
            xy = [x y]
        else
            xy = [a*(cos.(t) + 1/2pulse*cos.(pulse*t)) b*(sin.(t) + 1/2pulse*sin.(pulse*t))]
        end
    else
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

    return [P.M/2π * ones(N); thetas; 0; 0; 0]
end

end # module
