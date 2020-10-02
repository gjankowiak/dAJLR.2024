push!(LOAD_PATH, "src")

import Bend
import Plotting
import PyPlot
import Serialization
import Arpack


atol = 1e-8
rtol = 1e-13
max_iter = 3000
step_size = 1e-1

bigZ2(m, h) = h^3/2 - 9*h^2*m^2 + 18*h*m^4 - 7*m^6

solver_params = Bend.SolverParams(
                                  atol, rtol, max_iter, step_size,
                                  true, # adapt
                                  1e-5, # min step size
                                  1e-1, # max step size
                                  5e1, # step down threshold
                                  5e-4, # step up threshold
                                  1.4,  # step factor
                                 )

function iterate_filename(pattern)
    if occursin("#", pattern)
        for i in 1:1000
            s, s_P = string(replace(pattern, "#" => i), ".dat"), string(replace(pattern, "#" => i), "_P.dat")
            if (isfile(s) || isfile(s_P))
                continue
            end
            return s, s_P
        end
    else
        return string(pattern, ".dat"), string(pattern, "_P.dat")
    end
end

function follow_branch_any(P, param::Symbol, param_start, param_delta, n_samples; xinit="", prefix="result_#")
    Xs = []
    Ps = []

    if param_delta isa Number
        p = range(param_start, step=param_delta, length=n_samples)
    else
        p = Array{Float64,1}()
        e_last = param_start
        for (step, n) in param_delta
            if length(p) == 0
                p = collect(range(e_last, step=step, length=n))
            else
                p = [p; range(e_last+step, step=step, length=n)]
            end
            e_last = p[end]
        end
    end

    @show p

    if xinit isa String
        Xinit = Serialization.deserialize(xinit)[end]
    else
        Xinit = xinit
    end

    Xx = Base.copy(Xinit)

    f = Plotting.init(P)
    Plotting.plot(f, P, Xx)

    for (k,v) in enumerate(p)

        setproperty!(P, param, v)
        println(k, "/", length(p), " - ", string(param), ": ", getproperty(P, param))

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
                Plotting.plot(f, P, res.sol)
                push!(Xs, res.sol)
                push!(Ps, Bend.copy(P))
                break
            end
        end
    end

    s, s_P = iterate_filename(prefix)

    Serialization.serialize(s, Xs)
    Serialization.serialize(s_P, Ps)

    return Xs, Ps
end

function follow_branch(P, eps_delta, n_samples; eps_start=0.0, eps_direction=0.0, prefix="result_#", fn_xinit="")
    Xs = []
    Ps = []

    Z = bigZ2(P.beta_m, P.beta_h)
    eps_c = sqrt((P.beta_m^2 - P.beta_h/2)/P.mode_j^2)
    amp = sqrt(abs(P.mode_j^2*4*(2*P.beta_m^2 - P.beta_h)/Z))

    if eps_start > 0.0
        if eps_direction == 0.0
            throw("Please provide eps_direction")
        end
        epsilons = collect(eps_start .+ sign(eps_direction).*eps_delta*(0:n_samples-1))
    else
        if eps_delta isa Number
            epsilons = collect(eps_c .- sign(Z).*eps_delta*(1:n_samples))
        else
            p = Array{Float64,1}()
            e_last = eps_c
            for (step, n) in eps_delta
                p = [p; range(e_last-sign(Z)*step, step=-sign(Z)*step, length=n)]
                e_last = p[end]
            end
            epsilons = p
        end
    end

    @show epsilons
    P.epsilon = epsilons[1]

    @show Z
    println("j:", P.mode_j, ", eps: ", eps_c)
    println(", amplitude: ", amp)

    if fn_xinit == ""
        Xinit = Bend.initial_data(P, 1, 1, pulse=P.mode_j, pulse_amplitude=1e0*amp*sqrt(abs(eps_c - epsilons[1])), reverse_phase=false, only_rho=false)
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
        println(k, "/", length(epsilons), " - epsilon: ", P.epsilon)

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

    return Xs, Ps
end

function follow_branch_decoupled(P, eps_delta, n_samples; eps_start=0.0, eps_direction=0.0, prefix="result_#", fn_xinit="")
    Xs = []
    Ps = []

    if P.beta_m != 0 || P.beta_h >= 0
        throw("Parameter m is not zero or h is non negative")
    end
    eps_c = sqrt(-P.beta_h/2)
    amp = sqrt(8/P.beta_h^2)

    if eps_start > 0.0
        if eps_direction == 0.0
            throw("Please provide eps_direction")
        end
        epsilons = collect(eps_start .+ sign(eps_direction).*eps_delta*(0:n_samples-1))
    else
        if eps_delta isa Number
            epsilons = collect(eps_c .+ eps_delta*(1:n_samples))
        else
            p = Array{Float64,1}()
            e_last = eps_c
            for (step, n) in eps_delta
                p = [p; range(e_last, step=step, length=n)]
                @show p
                e_last = p[end]
            end
            epsilons = p
        end
    end

    @show epsilons


    println("j:", P.mode_j, ", eps: ", eps_c)
    println(", amplitude: ", amp)

    if fn_xinit == ""
        Xinit = Bend.initial_data(P, 1, 1, pulse=P.mode_j, pulse_amplitude=amp*sqrt(abs(eps_c - epsilons[1])), reverse_phase=false, only_rho=false)
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
                Plotting.plot(f, P, res.sol)
                break
            end
        end
    end

    s, s_P = iterate_filename(prefix)

    Serialization.serialize(s, Xs)
    Serialization.serialize(s_P, Ps)

    return Xs, Ps
end

function find_critical_point(P, eps_delta)
    Xs = []
    Ps = []

    eps_c = sqrt((P.beta_m^2 - P.beta_h/2)/P.mode_j^2)
    Z = bigZ2(P.beta_m, P.beta_h)
    amp = sqrt(abs(P.mode_j^2*4*(2*P.beta_m^2 - P.beta_h)/Z))

    println("2m^2 - h = ", 2*P.beta_m^2 - P.beta_h)

    print("j:", P.mode_j, ", eps: ", eps_c)
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

#include("run_bend_old_params.jl")

# include("params/check_case_0_stabilty.jl")

# include("params/continue_case_1.1.jl")
#
# include("params/continue_case_1.1_m-1h-1.jl")

include("params/case_0_m=1.jl")

# include("params/case_0_m=-1.jl")
