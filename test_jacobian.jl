push!(LOAD_PATH, "src")

import Bend
import PyPlot
import LinearAlgebra
const LA = LinearAlgebra

N = 600
Δs = 2π/N
M = 2π

epsilon = 1e-1
potential_range = 0
rho_max = 3*M/2π

center_rho = false

atol = 1e-3
rtol = 1e-13
max_iter = 50000
step_size = 1e-3

beta_0 = 1
beta_rho0  = M/2π
beta_m     = -2
beta_h     = -2
beta_k     = 20
beta_j     = 4

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
                                  1e-3, # step up threshold
                                  1.2,  # step factor
                                 )

test = false
if test
    c_r = Bend.compute_residual_test
    c_J = Bend.assemble_inner_system_test
else
    c_r = Bend.compute_residual
    c_J = Bend.assemble_inner_system
end

Xinit = Bend.initial_data(P, 1, 1, pulse=4, reverse_phase=true)
Xx = copy(Xinit)

h_i = 10. .^range(-8,0,length=10)
e_i = zeros(size(h_i))
ebeta = zeros(size(h_i))


matrices = Bend.assemble_fd_matrices(P)

mask_in = [ones(P.N); ones(P.N-1); ones(3)]
mask_out = [ones(P.N); ones(P.N-1); ones(3)]

r0 = c_r(P, matrices, Xx)
d = rand(size(r0)...) .- 0.5
J0 = c_J(P, matrices, Xx)

xbeta = 1.4
r0beta = Bend.compute_beta_prime(P, [xbeta])
J0beta = Bend.compute_beta_second(P, [xbeta])

for (i,h) in enumerate(h_i)
    rh = c_r(P, matrices, Xx + h*(mask_in.*d))
    e_i[i] = LA.norm(mask_out.*(r0 + J0*(h*(mask_in.*d)) - rh))

    rhbeta = Bend.compute_beta_prime(P, [xbeta + h])
    ebeta[i] = LA.norm(r0beta + J0beta*h - rhbeta)
end

slope = (log(e_i[3]) - log(e_i[end-3])) / (log(1/h_i[3]) - log(1/h_i[end-3]))
slope2 = (log(ebeta[1]) - log(ebeta[end])) / (log(1/h_i[1]) - log(1/h_i[end]))

PyPlot.subplot(121)
PyPlot.loglog(1 ./h_i, e_i)
PyPlot.title(slope)
PyPlot.grid()
PyPlot.subplot(122)
PyPlot.loglog(1 ./h_i, ebeta)
PyPlot.title(slope2)
PyPlot.grid()
PyPlot.show()
