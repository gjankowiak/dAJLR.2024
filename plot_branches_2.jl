push!(LOAD_PATH, "src")

import Bend
import PyPlot
import Serialization
import Interpolations

color = ["#1f77b4", "#ff7f0e"]

case = 1
function plotall(case::Int)
    if case == 1
        fn_Xs_b1 = "results_case_1/branch_1.dat"
        fn_Ps_b1 = "results_case_1/branch_1_P.dat"
        fn_Xs_b2 = "results_case_1/branch_2.dat"
        fn_Ps_b2 = "results_case_1/branch_2_P.dat"
        eps_c = [0.353553390593273; 0.23570226039551584]
        title = "β = 1 + (ρ-1) + 1/2 (ρ-1)²"
    elseif case == 2
        fn_Xs_b1 = "results_case_2/branch_1.dat"
        fn_Ps_b1 = "results_case_2/branch_1_P.dat"
        fn_Xs_b2 = "results_case_2/branch_2.dat"
        fn_Ps_b2 = "results_case_2/branch_2_P.dat"
        eps_c = [1.118033988749895; 0.7453559924999299]
        title = "β = 1 - 2 (ρ-1) - (ρ-1)² + 20 (ρ-1)^4"
    elseif case == 3
        fn_Xs_b1 = "results_case_3/branch_1.dat"
        fn_Ps_b1 = "results_case_3/branch_1_P.dat"
        fn_Xs_b2 = "results_case_3/branch_2.dat"
        fn_Ps_b2 = "results_case_3/branch_2_P.dat"
        eps_c = [0.6123724356957945; 0.408248290463863]
        title = "β = 1 - 1 (ρ-1) - 1/2 (ρ-1)² + 20 (ρ-1)^4"
    elseif case == 4
        fn_Xs_b1 = "results_case_4/branch_1.dat"
        fn_Ps_b1 = "results_case_4/branch_1_P.dat"
        fn_Xs_b2 = "results_case_4/branch_2.dat"
        fn_Ps_b2 = "results_case_4/branch_2_P.dat"
        eps_c = [0.3535533905932738; 0.23570226039551584]
        title = "β = 1 - 1 (ρ-1) + 1/2 (ρ-1)²"
    elseif case == 5
        fn_Xs_b1 = "results_case_5/branch_1.dat"
        fn_Ps_b1 = "results_case_5/branch_1_P.dat"
        fn_Xs_b2 = "results_case_5/branch_1.dat"
        fn_Ps_b2 = "results_case_5/branch_1_P.dat"
        eps_c = [0.13635890143294643]
        title = "β = 1 - 1.3 (ρ-1) + 1/2 (ρ-1)² + lower bound on ρ"
    else
        error("case not done")
    end

    # branch 1
    Xs_b1 = Serialization.deserialize(fn_Xs_b1)
    Ps_b1 = Serialization.deserialize(fn_Ps_b1)

    epsilons_b1 = map(x -> x.epsilon, Ps_b1)
    mass_energy_b1 = zeros(size(epsilons_b1))
    bending_energy_b1 = zeros(size(epsilons_b1))
    min_rho_b1 = zeros(size(epsilons_b1))
    max_rho_b1 = zeros(size(epsilons_b1))


    # branch 2
    Xs_b2 = Serialization.deserialize(fn_Xs_b2)
    Ps_b2 = Serialization.deserialize(fn_Ps_b2)
    epsilons_b2 = map(x -> x.epsilon, Ps_b2)
    mass_energy_b2 = zeros(size(epsilons_b2))
    bending_energy_b2 = zeros(size(epsilons_b2))
    min_rho_b2 = zeros(size(epsilons_b2))
    max_rho_b2 = zeros(size(epsilons_b2))

    matrices = Bend.assemble_fd_matrices(Ps_b1[1])

    P = Ps_b1[1]
    N = P.N
    Δs = P.Δs

    t = collect(range(0, 2π, length=N+1))[1:N]

    xy1 = zeros(N+1, 2)
    xy2 = zeros(N+1, 2)
    xy1[2,:] = [Δs 0]

    for i in 1:length(Xs_b1)
        P = Ps_b1[i]
        (mass_energy_b1[i], bending_energy_b1[i]) = Bend.compute_energy_split(P, matrices, Xs_b1[i])

        c = Bend.X2candidate(P, Xs_b1[i])
        min_rho_b1[i], max_rho_b1[i] = extrema(c.ρ)
    end

    for i in 1:length(Xs_b2)
        P = Ps_b2[i]
        (mass_energy_b2[i], bending_energy_b2[i]) = Bend.compute_energy_split(P, matrices, Xs_b2[i])

        c = Bend.X2candidate(P, Xs_b2[i])
        min_rho_b2[i], max_rho_b2[i] = extrema(c.ρ)
    end

    fig = PyPlot.figure(dpi=50)
    ax1 = PyPlot.subplot(241)
    ax1.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
    ax1.plot(epsilons_b1, mass_energy_b1.+bending_energy_b1, label="branch 1", ".-")
    ax1.plot(epsilons_b2, mass_energy_b2.+bending_energy_b2, label="branch 2", ".-")

    b1_intercept = (((mass_energy_b1[1]+bending_energy_b1[1])*epsilons_b1[2] -
                     (mass_energy_b1[2]+bending_energy_b1[2])*epsilons_b1[1]) /
                    (epsilons_b1[2] - epsilons_b1[1]))

    b2_intercept = (((mass_energy_b2[1]+bending_energy_b2[1])*epsilons_b2[2] -
                     (mass_energy_b2[2]+bending_energy_b2[2])*epsilons_b2[1]) /
                    (epsilons_b2[2] - epsilons_b2[1]))

    ax1.plot([0; epsilons_b1[1]], [b1_intercept; mass_energy_b1[1]+bending_energy_b1[1]], color="black")
    ax1.plot([0; epsilons_b2[1]], [b2_intercept; mass_energy_b2[1]+bending_energy_b2[1]], color="black")

    ax1.axhline(π, label="circle", color="black")
    ax1.axvline(0, color="black", lw=0.5)
    for eps in eps_c
        ax1.axvline(eps, color="black", lw=0.5)
    end
    ax1.set_title("Energy")
    ax1.legend()

    ax2 = PyPlot.subplot(245, sharex=ax1)
    ax2.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
    ax2.plot(epsilons_b1, mass_energy_b1, ".-", color=color[1])
    ax2.plot(epsilons_b1, bending_energy_b1, ".:",  color=color[1])
    ax2.plot(epsilons_b2, mass_energy_b2, ".-",  color=color[2])
    ax2.plot(epsilons_b2, bending_energy_b2, ".:",  color=color[2])
    ax2.axhline(0, color="black", lw=0.5)
    ax2.axvline(0, color="black", lw=0.5)
    ax2.set_title("Energy components")

    ax3 = PyPlot.subplot(142)
    ax3.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
    ax3.plot(epsilons_b1, min_rho_b1, label="branch 1", ".-", color=color[1])
    ax3.plot(epsilons_b1, max_rho_b1, label="branch 1", ".-", color=color[1])
    ax3.plot(epsilons_b2, min_rho_b2, label="branch 1", ".-", color=color[2])
    ax3.plot(epsilons_b2, max_rho_b2, label="branch 1", ".-", color=color[2])
    ax3.axhline(1.0, label="circle", color="black")
    for eps in eps_c
        ax3.axvline(eps, color="black", lw=0.5)
    end
    ax3.set_title("Min/max rho")
    ax3.legend()

    ax4 = PyPlot.subplot(243)
    ax5 = PyPlot.subplot(247)

    ax6 = PyPlot.subplot(4, 4, 4)
    ax7 = PyPlot.subplot(4, 4, 8)
    ax8 = PyPlot.subplot(4, 4, 12)
    ax9 = PyPlot.subplot(4, 4, 16)

    findnearest(A::Vector{Float64},t::Float64) = findmin(abs.(A.-t))[2]

    function plot_curve(i1::Int64, i2::Int64; init::Bool=false)
        X1 = Xs_b1[i1]
        X2 = Xs_b2[i2]
        ax = [ax4, ax5]
        ax_rho = [ax6, ax8]
        ax_theta = [ax7, ax9]
        xy_pos = [xy1, xy2]
        for (k,X) in enumerate((X1, X2))
            xy = xy_pos[k]
            c = Bend.X2candidate(P, X)

            for i in 3:N+1
                xy[i,:] = xy[i-1,:] + Δs*[cos(c.θ[i-2]); sin(c.θ[i-2])]
            end

            bc = sum(xy, dims=1)/N
            xy .= xy .- bc

            θ_dot = Bend.compute_centered_fd_θ(P, matrices, c.θ)

            if init
                ax[k].plot(xy[:,1], xy[:,2], lw=0.1, color=color[k])
                ax[k].scatter(xy[:,1], xy[:,2], s=3*c.ρ, color=color[k])
                ax[k].set_aspect("equal")
                ax[k].set_xlim(-1.5, 1.5)
                ax[k].set_ylim(-1.5, 1.5)
                ax[k].set_title("2D visualisation")

                ax_rho[k].plot(t, c.ρ, color=color[k])
                ax_rho[k].axhline(0, color="black", lw=0.5)
                ax_rho[k].set_title("ρ")

                ax_theta[k].plot(t, θ_dot, color=color[k])
                ax_theta[k].set_title("d/ds θ")
                ax_theta[k].set_ylim(0, 10)
            else
                ax[k].lines[1].set_data(xy[:,1], xy[:,2])
                ax[k].collections[1].set_offsets(xy)
                ax[k].collections[1].set_sizes(3*c.ρ)

                ax_rho[k].lines[1].set_data(t, c.ρ)
                ax_theta[k].lines[1].set_data(t, θ_dot)

                ax_rho[k].relim()
                ax_rho[k].autoscale_view(true,true,true)
            end
            PyPlot.draw()
        end
    end

    plot_curve(1, 1; init=true)

    function mousemove(event)
        if event.inaxes in [ax1, ax2, ax3]
            for ax in [ax1, ax2, ax3]
                d = ax.lines[1].get_data()
                ax.lines[1].set_data([event.xdata], d[2])
            end
        end
    end

    function mouseup(event)
        if event.inaxes in [ax1, ax2, ax3]
            idx1 = findnearest(epsilons_b1, event.xdata)
            idx2 = findnearest(epsilons_b2, event.xdata)
            plot_curve(idx1, idx2)
        end
    end

    fig.canvas.mpl_connect("motion_notify_event", mousemove)
    fig.canvas.mpl_connect("button_release_event", mouseup)

    fig.suptitle(title)

    PyPlot.show()
end

# plotall(1)
# plotall(2)
# plotall(3)
plotall(4)
# plotall(5)
