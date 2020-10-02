push!(LOAD_PATH, "src")

import Bend
import PyPlot
import Serialization
import Interpolations

color = ["#1f77b4", "#ff7f0e"]

eps_critical(m, h, j) = sqrt((m^2 - 0.5*h)/j^2)
bigZ2(m, h) = 0.5*(h^3 - 18*h^2*m^2 + 36*h*m^4 - 14*m^6)
amp2(m, h, j) = sqrt(abs(j^2*4*(2*m^2 - h)/bigZ2(m, h)))

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
        eps_c = [0.13635890143294643, 0]
        title = "β = 1 - 1.3 (ρ-1) + 1/2 (ρ-1)² + lower bound on ρ"
    elseif case == 7
        fn_Xs_b1 = "results_case_7/branch_1.dat"
        fn_Ps_b1 = "results_case_7/branch_1_P.dat"
        fn_Xs_b2 = "results_case_7/branch_2.dat"
        fn_Ps_b2 = "results_case_7/branch_2_P.dat"
        # fn_Xs_b2 = "results_case_7/branch_1_rerun.dat"
        # fn_Ps_b2 = "results_case_7/branch_1_rerun_P.dat"
        eps_c = [0.6123724356957945; 0.408248290463863]
        title = "β = 1 - (ρ-1) - 1/2 (ρ-1)²"
    elseif case == 8
        fn_Xs_b1 = "results_case_8/branch_1_new.dat"
        fn_Ps_b1 = "results_case_8/branch_1_new_P.dat"
        fn_Xs_b2 = "results_case_8/branch_2_new.dat"
        fn_Ps_b2 = "results_case_8/branch_2_new_P.dat"
        eps_c = [0.726291952; 0.48419463]
        title = "β = 1 + 1.9(ρ-1) + 3/2 (ρ-1)²"
    elseif case == 9
        fn_Xs_b1 = "results_decoupled/branch_0.dat"
        fn_Ps_b1 = "results_decoupled/branch_0_P.dat"
        fn_Xs_b2 = "results_decoupled/branch_1.dat"
        fn_Ps_b2 = "results_decoupled/branch_1_P.dat"
        # fn_Xs_b2 = "results_decoupled/branch_2.dat"
        # fn_Ps_b2 = "results_decoupled/branch_2_P.dat"
        eps_c = [0.7071067811865476, 0.3535533905932738]
        # eps_c = [0.3535533905932738; 0.23570226039551584]
        title = "β = 1 - 1/2 (ρ-1)²"
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
    ax1.set_ylabel("Energy")
    ax1.legend()

    ax2 = PyPlot.subplot(245, sharex=ax1)
    ax2.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
    ax2.plot(epsilons_b1, mass_energy_b1, color=color[1], label="\$\\frac{1}{2}\\int \\dot{\\rho}^2\$")
    ax2.plot(epsilons_b1, bending_energy_b1, ":",  color=color[1], label="\$\\frac{1}{2}\\int \\beta(\\rho) \\dot{\\theta}^2\$")
    ax2.plot(epsilons_b2, mass_energy_b2, color=color[2])
    ax2.plot(epsilons_b2, bending_energy_b2, ":",  color=color[2])
    ax2.axhline(0, color="black", lw=0.5)
    ax2.axvline(0, color="black", lw=0.5)
    ax2.legend()
    ax2.set_ylabel("Energy components")

    ax3 = PyPlot.subplot(142)
    ax3.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
    ax3.plot(epsilons_b1, min_rho_b1, label="branch 1", ".-", color=color[1])
    ax3.plot(epsilons_b1, max_rho_b1, label="branch 1", ".-", color=color[1])
    ax3.plot(epsilons_b2, min_rho_b2, label="branch 1", ".-", color=color[2])
    ax3.plot(epsilons_b2, max_rho_b2, label="branch 1", ".-", color=color[2])

    Z = bigZ2(P.beta_a1, P.beta_a2)
    @show P.beta_a1
    @show P.beta_a2
    @show Z

    amp_b1 = sqrt.(ifelse.((sign(Z).*(eps_c[1] .- epsilons_b1) .> 0), sign(Z).*(eps_c[1] .- epsilons_b1), NaN)).*amp2(P.beta_a1, P.beta_a2, 2)
    amp_b2 = sqrt.(ifelse.((sign(Z).*(eps_c[2] .- epsilons_b2) .> 0), sign(Z).*(eps_c[2] .- epsilons_b2), NaN)).*amp2(P.beta_a1, P.beta_a2, 3)

    println(size(amp_b1))
    println(size(epsilons_b1))
    for s in [-1, 1]
        ax3.plot(epsilons_b1, 1 .+ s*amp_b1, ls="dotted", color="black", lw=0.5)
        ax3.plot(epsilons_b2, 1 .+ s*amp_b2, ls="dotted", color="black", lw=0.5)
    end

    ax3.axhline(1.0, label="circle", color="black")
    for j in [2, 3]
        ax3.axvline(eps_critical(P.beta_a1, P.beta_a2, j), color="black", lw=0.5)
    end
    ax3.set_ylabel("Min/max rho")
    ax3.legend()

    ax4 = PyPlot.subplot(243)
    ax5 = PyPlot.subplot(247)

    ax6 = PyPlot.subplot(4, 4, 4)
    ax6bis = ax6.twinx()
    ax7 = PyPlot.subplot(4, 4, 8)
    ax8 = PyPlot.subplot(4, 4, 12)
    ax8bis = ax8.twinx()
    ax9 = PyPlot.subplot(4, 4, 16)

    findnearest(A::Vector{Float64},t::Float64) = findmin(abs.(A.-t))[2]

    function plot_curve(i1::Int64, i2::Int64; init::Bool=false)
        X1 = Xs_b1[i1]
        X2 = Xs_b2[i2]
        ax = [ax4, ax5]
        ax_rho = [ax6, ax6bis, ax8, ax8bis]
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

            beta = Bend.compute_beta(P, c.ρ)

            θ_dot = Bend.compute_centered_fd_θ(P, matrices, c.θ)

            if init
                ax[k].plot(xy[:,1], xy[:,2], lw=0.1, color=color[k])
                ax[k].scatter(xy[:,1], xy[:,2], s=3*c.ρ, color=color[k])
                ax[k].set_aspect("equal")
                ax[k].set_xlim(-1.5, 1.5)
                ax[k].set_ylim(-1.5, 1.5)
                # ax[k].set_title("2D visualisation")

                ax_rho[2*k-1].plot(t, c.ρ, color=color[k])
                ax_rho[2*k-1].axhline(0, color="black", lw=0.5)
                ax_rho[2*k].plot(t, beta, ls="dotted", color=color[k])
                ax_rho[2*k-1].set_ylabel("ρ")
                ax_rho[2*k].set_ylabel("β(ρ)")

                ax_theta[k].plot(t, θ_dot, color=color[k])
                ax_theta[k].set_ylabel("d/ds θ")
                ax_theta[k].set_ylim(0, 10)
            else
                ax[k].lines[1].set_data(xy[:,1], xy[:,2])
                ax[k].collections[1].set_offsets(xy)
                ax[k].collections[1].set_sizes(3*c.ρ)

                ax_rho[2*k-1].lines[1].set_data(t, c.ρ)
                ax_rho[2*k].lines[1].set_data(t, beta)
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
end

# plotall(1)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(2)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(3)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(4)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(5)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(7)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

# plotall(8)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

plotall(9)
mng = PyPlot.get_current_fig_manager()
mng.window.showMaximized()
sleep(1)
PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

PyPlot.show()
