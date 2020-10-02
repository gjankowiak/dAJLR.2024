push!(LOAD_PATH, "src")

import Bend
# import Plotting
# import PyPlot
import Serialization

fn_X = ARGS[1]
fn_P = replace(fn_X, ".dat" => "_P.dat")

fn_X_new = replace(fn_X, ".dat" => "_rerun.dat")
fn_P_new = replace(fn_P, "_P.dat" => "_rerun_P.dat")

Xs_old = Serialization.deserialize(fn_X)
Ps_old = Serialization.deserialize(fn_P)

Xs = []
Ps = []

atol = 1e-8
rtol = 1e-13
max_iter = 3000
step_size = 1e-2

solver_params = Bend.SolverParams(
                                  atol, rtol, max_iter, step_size,
                                  true, # adapt
                                  1e-5, # min step size
                                  1e-1, # max step size
                                  5e1, # step down threshold
                                  5e-4, # step up threshold
                                  1.4,  # step factor
                                 )

Xx = zeros(size(Xs_old[1]))

for (i,P) in enumerate(Ps_old)

    Xx .= Xs_old[i]

    # change parameters
    # P.beta_a1 = 2e-1
    # P.beta_a2 = -1.8
    # P.beta_a4 = 2

    # f = Plotting.init(P)
    # Plotting.plot(f, P, Xx)

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
            println(res.residual_norm_i[n])
            # Plotting.plot(f, P, res.sol)
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

Serialization.serialize(fn_X_new, Xs)
Serialization.serialize(fn_P_new, Ps)

# for r in Xs
    # Plotting.plot(P, r)
# end

# display(X)
# println()
# println(err)
