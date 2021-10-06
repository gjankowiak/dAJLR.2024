push!(LOAD_PATH, "src")

import Bend
import PyPlot
import Serialization
import Interpolations

import DelimitedFiles

color = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]

eps_critical(m, h, j) = sqrt((m^2 - 0.5*h)/j^2)
bigZ2(m, h) = 0.5*(h^3 - 18*h^2*m^2 + 36*h*m^4 - 14*m^6)
amp2(m, h, j) = sqrt(abs(j^2*4*(2*m^2 - h)/bigZ2(m, h)))

case = 1
function plotall(case::Int; only_branch=0, dot_size_factor::Float64=3.0)
    if case == 1
        fns_Xs = ["results_case_1/branch_1.dat", "results_case_1/branch_2.dat"]
        fns_Ps = ["results_case_1/branch_1_P.dat", "results_case_1/branch_2_P.dat"]
        title = "β = 1 + (ρ-1) + 1/2 (ρ-1)²"
    elseif case == 2
        fns_Xs = ["results_case_2/branch_1.dat", "results_case_2/branch_2.dat"]
        fns_Ps = ["results_case_2/branch_1_P.dat", "results_case_2/branch_2_P.dat"]
        title = "β = 1 - 2 (ρ-1) - (ρ-1)² + 20 (ρ-1)^4"
    elseif case == 3
        fns_Xs = ["results_case_3/branch_1.dat", "results_case_3/branch_2.dat"]
        fns_Ps = ["results_case_3/branch_1_P.dat", "results_case_3/branch_2_P.dat"]
        title = "β = 1 - 1 (ρ-1) - 1/2 (ρ-1)² + 20 (ρ-1)^4"
    elseif case == 4
        fns_Xs = ["results_case_4/branch_1.dat", "results_case_4/branch_2.dat"]
        fns_Ps = ["results_case_4/branch_1_P.dat", "results_case_4/branch_2_P.dat"]
        title = "β = 1 - 1 (ρ-1) + 1/2 (ρ-1)²"
    elseif case == 5
        fns_Xs = ["results_case_5/branch_1.dat"]
        fns_Ps = ["results_case_5/branch_1_P.dat"]
        title = "β = 1 - 1.3 (ρ-1) + 1/2 (ρ-1)² + lower bound on ρ"
    elseif case == 7
        fns_Xs = ["results_case_7/branch_1.dat", "results_case_7/branch_2.dat"]
        fns_Ps = ["results_case_7/branch_1_P.dat", "results_case_7/branch_2_P.dat"]
        # fn_Xs_b2 = "results_case_7/branch_1_rerun.dat"
        # fn_Ps_b2 = "results_case_7/branch_1_rerun_P.dat"
        title = "β = 1 - (ρ-1) - 1/2 (ρ-1)²"
    elseif case == 8
        fns_Xs = ["results_case_8/branch_1.dat", "results_case_8/branch_2.dat"]
        fns_Ps = ["results_case_8/branch_1_P.dat", "results_case_8/branch_2_P.dat"]
        title = "β = 1 + 1.9(ρ-1) + 3/2 (ρ-1)²"
    elseif case == 9
        fns_Xs = ["results_decoupled/branch_0.dat", "results_decoupled/branch_1.dat", "results_decoupled/branch_2.dat"]
        fns_Ps = ["results_decoupled/branch_0_P.dat", "results_decoupled/branch_1_P.dat", "results_decoupled/branch_2_P.dat"]
        # fns_Xs = ["results_decoupled/branch_0_rerun.dat",
                  # "results_decoupled/branch_0.dat"]
                  # #"results_decoupled/branch_2_rerun.dat"]
        # fns_Ps = ["results_decoupled/branch_0_rerun_P.dat",
                  # "results_decoupled/branch_0_P.dat"]
                  # #"results_decoupled/branch_2_rerun_P.dat"]
        title = "β = 1 - 1/2 (ρ-1)²"
    elseif case == 10
        fns_Xs = ["Pm2m2/branch_1.dat", "Pm2m2/branch_2.dat"]
        fns_Ps = ["Pm2m2/branch_1_P.dat", "Pm2m2/branch_2_P.dat"]
        title = "β = 1 - 2 (ρ-1) - 2/2 (ρ-1)²"
    elseif case == 100
        fns_Xs = ["results_check_stability/check_stability_sub_stable_1_new.dat", "results_check_stability/check_stability_super_unstable_1_new.dat"]
        fns_Ps = ["results_check_stability/check_stability_sub_stable_1_new_P.dat", "results_check_stability/check_stability_super_unstable_1_new_P.dat"]
        title = "Case 0 Stability"
    elseif case == 101
        fns_Xs = ["results_continue_case_1.1/continue_case_1.1_sub_c_1.dat", "results_continue_case_1.1/continue_case_1.1_super_c_1.dat", "results_continue_case_1.1/continue_case_1.1_super_c2_1.dat"]
        fns_Ps = ["results_continue_case_1.1/continue_case_1.1_sub_c_1_P.dat", "results_continue_case_1.1/continue_case_1.1_super_c_1_P.dat", "results_continue_case_1.1/continue_case_1.1_super_c2_1_P.dat"]
        title = "Case 1.1 Stability"
    elseif case == 102
        fns_Xs = ["case0_m=1/stable_j=2_1.dat", "case0_m=1/unstable_low_j=2_1.dat", "case0_m=1/unstable_high_j=2_1.dat"]
        fns_Ps = ["case0_m=1/stable_j=2_1_P.dat", "case0_m=1/unstable_low_j=2_1_P.dat", "case0_m=1/unstable_high_j=2_1_P.dat"]
        title = "Comparison m=1, j=2"
    elseif case == 103
        fns_Xs = ["case0_m=1/stable_j=3_1.dat", "case0_m=1/unstable_low_j=3_1.dat", "case0_m=1/unstable_high_j=3_1.dat"]
        fns_Ps = ["case0_m=1/stable_j=3_1_P.dat", "case0_m=1/unstable_low_j=3_1_P.dat", "case0_m=1/unstable_high_j=3_1_P.dat"]
        title = "Comparison m=1, j=3"
    elseif case == 104
        fns_Xs = ["case11_m=1/stable_m=crit_pre_eps_1.dat", "case11_m=1/stable_m=crit_pre_eps_c0_j2_1.dat", "case11_m=1/stable_m=crit_pre_eps_c0_j3_1.dat", "case11_m=1/stable_m=crit_post_eps_1.dat"]
        fns_Ps = ["case11_m=1/stable_m=crit_pre_eps_1_P.dat", "case11_m=1/stable_m=crit_pre_eps_c0_j2_1_P.dat", "case11_m=1/stable_m=crit_pre_eps_c0_j3_1_P.dat", "case11_m=1/stable_m=crit_post_eps_1_P.dat"]
        title = "Case 11 Comparison m=1, j=3"
    elseif case == 201 # unstable high Case 0
        fns_Xs = ["case0_m=1/unstable_high_j=2_1_new.dat", "case0_m=1/unstable_high_j=3_1_new.dat", "case0_m=1/unstable_high_j=1_1.dat"]
        fns_Ps = ["case0_m=1/unstable_high_j=2_1_new_P.dat", "case0_m=1/unstable_high_j=3_1_new_P.dat", "case0_m=1/unstable_high_j=1_1_P.dat"]
        title = "Parameter set (i)"
    elseif case == 202 # stable Case 0
        fns_Xs = ["case0_m=1/stable_j=2_1.dat", "case0_m=1/stable_j=3_1.dat", "case0_m=1/stable_j=1_1.dat"]
        fns_Ps = ["case0_m=1/stable_j=2_1_P.dat", "case0_m=1/stable_j=3_1_P.dat", "case0_m=1/stable_j=1_1_P.dat"]
        title = "Parameter set (ii)"
    elseif case == 203 # unstable low Case 0
        fns_Xs = ["case0_m=1/unstable_low_j=2_1.dat", "case0_m=1/unstable_low_j=3_1.dat", "case0_m=1/unstable_low_j=1_1.dat"]
        fns_Ps = ["case0_m=1/unstable_low_j=2_1_P.dat", "case0_m=1/unstable_low_j=3_1_P.dat", "case0_m=1/unstable_low_j=1_1_P.dat"]
        title = "Parameter set (iii)"
    elseif case == 204 # stable Case 1.1
        fns_Xs = ["case11_m=1/stable_m=1_eps_1_new.dat",   "case11_m=1/stable_m=1_eps_c0_j2_1.dat",   "case11_m=1/stable_m=1_eps_c0_j3_1.dat"]
        fns_Ps = ["case11_m=1/stable_m=1_eps_1_new_P.dat", "case11_m=1/stable_m=1_eps_c0_j2_1_P.dat", "case11_m=1/stable_m=1_eps_c0_j3_1_P.dat"]
        title = "Parameter set (iv)"
    elseif case == 205 # unstable Case 1.1
        fns_Xs = ["case11_m=1/unstable_m=1_eps_1.dat", "case11_m=1/unstable_m=1_eps_c0_j2_1.dat", "case11_m=1/unstable_m=1_eps_c0_j3_1.dat"]
        fns_Ps = ["case11_m=1/unstable_m=1_eps_1_P.dat", "case11_m=1/unstable_m=1_eps_c0_j2_1_P.dat", "case11_m=1/unstable_m=1_eps_c0_j3_1_P.dat"]
        title = "Parameter set (v)"
    else
        error("case not done")
    end

    n_branches = length(fns_Xs)

    l_Xs = Dict{Int64,Array{Vector{Float64},1}}()
    l_Ps = Dict{Int64,Array{Bend.ParamsUnion,1}}()

    l_epsilons = Dict{Int64, Vector{Float64}}()
    l_mass_energy = Dict{Int64, Vector{Float64}}()
    l_bending_energy = Dict{Int64, Vector{Float64}}()
    l_min_rho = Dict{Int64, Vector{Float64}}()
    l_max_rho = Dict{Int64, Vector{Float64}}()
    l_min_beta = Dict{Int64, Vector{Float64}}()
    l_max_beta = Dict{Int64, Vector{Float64}}()

    l_stab = Dict{Int64, Vector{Float64}}()

    l_eps_c = zeros(n_branches)

    local matrices, xy, t, P

    for b in 1:n_branches
        l_Xs[b] = Serialization.deserialize(fns_Xs[b])
        l_Ps[b] = Serialization.deserialize(fns_Ps[b])

        # stability
        fn_stab = replace(fns_Xs[b], (".dat" => "_flowed_errors.dat"))
        if isfile(fn_stab)
            stab = Serialization.deserialize(fn_stab)
            l_stab[b] = abs.(stab[:,3] - stab[:,4])
        else
            l_stab[b] = ones(size(l_Ps[b], 1))
        end

        if b == 1
            matrices = Bend.assemble_fd_matrices(l_Ps[1][1])
            xy = zeros(l_Ps[1][1].N+1, 2)
            t = collect(range(0, 2π, length=l_Ps[1][1].N+1))[1:l_Ps[1][1].N]
        end

        l_epsilons[b] = map(x -> x.epsilon, l_Ps[b])

        n_samples = size(l_epsilons[b])

        l_mass_energy[b]    = zeros(n_samples)
        l_bending_energy[b] = zeros(n_samples)
        l_min_rho[b]        = zeros(n_samples)
        l_max_rho[b]        = zeros(n_samples)
        l_min_beta[b]        = zeros(n_samples)
        l_max_beta[b]        = zeros(n_samples)

        P = l_Ps[b][1]
        l_eps_c[b] = sqrt((P.beta_a1^2 - P.beta_a2/2)/P.mode_j^2)

        for i in 1:n_samples[1]
            P = l_Ps[b][i]
            (l_mass_energy[b][i], l_bending_energy[b][i]) = Bend.compute_energy_split(P, matrices, l_Xs[b][i])

            c = Bend.X2candidate(l_Ps[b][i], l_Xs[b][i])
            # FIXME
            l_min_rho[b][i], l_max_rho[b][i] = extrema(c.ρ)
            l_min_beta[b][i], l_max_beta[b][i] = extrema(Bend.compute_beta(l_Ps[b][i], c.ρ))
        end

        data = [l_epsilons[b] l_mass_energy[b] .+ l_bending_energy[b] l_mass_energy[b] l_bending_energy[b] l_min_rho[b] l_max_rho[b] l_min_beta[b] l_max_beta[b]]
        DelimitedFiles.writedlm("Data_case_$(case)_branch_$(b).csv", data)
    end


    fig = PyPlot.figure(dpi=50)

    ax1 = PyPlot.subplot(222)
    ax1.axvline(l_epsilons[1][1], ls="dotted", color="black", lw=0.5)
    ax1.axhline(π, label="circle", color="black")
    ax1.axvline(0, color="black", lw=0.5)
    ax1.set_ylabel("Energy")

    for b in 1:n_branches
        ax1.plot(l_epsilons[b], l_mass_energy[b].+l_bending_energy[b], label=string("branch ", b), ".-")

        b_intercept = (((l_mass_energy[b][1]+l_bending_energy[b][1])*l_epsilons[b][2] -
                         (l_mass_energy[b][2]+l_bending_energy[b][2])*l_epsilons[b][1]) /
                        (l_epsilons[b][2] - l_epsilons[b][1]))

        ax1.plot([0; l_epsilons[b][1]], [b_intercept; l_mass_energy[b][1]+l_bending_energy[b][1]], color="black")

        ax1.axvline(l_eps_c[b], color="black", lw=0.5)
    end
    ax1.legend()

    ax2 = PyPlot.subplot(224, sharex=ax1)
    ax2.axvline(l_epsilons[1][1], ls="dotted", color="black", lw=0.5)
    ax2.axhline(0, color="black", lw=0.5)
    ax2.axvline(0, color="black", lw=0.5)
    ax2.set_ylabel("Energy components")

    for b in 1:n_branches
        ax2.axvline(l_eps_c[b], color="black", lw=0.5)
        label1 = if b == 1 "\$\\frac{\\varepsilon^2}{2}\\int \\dot{\\rho}^2\$" else "" end
        label2 = if b == 1 "\$\\frac{1}{2}\\int \\beta(\\rho) \\dot{\\theta}^2\$" else "" end
        ax2.plot(l_epsilons[b], l_mass_energy[b], color=color[b], label=label1)
        ax2.plot(l_epsilons[b], l_bending_energy[b], ":",  color=color[b], label=label2)
    end

    ax2.legend()

    ax3 = PyPlot.subplot(121)


    ax3.axvline(l_epsilons[1][1], ls="dotted", color="black", lw=0.5)
    ax3.axhline(1.0, label="circle", color="black")
    ax3.set_ylabel("Min/max rho")
    for b in 1:n_branches
        ax3.scatter(l_epsilons[b], l_min_rho[b], s=100*(l_stab[b] .< 1e-7), color=color[b])
        ax3.scatter(l_epsilons[b], l_max_rho[b], s=100*(l_stab[b] .< 1e-7), color=color[b])
        ax3.axvline(l_eps_c[b], color="black", lw=0.5)
        ax3.plot(l_epsilons[b], l_min_rho[b], label=string("branch ", b), ".-", color=color[b])
        ax3.plot(l_epsilons[b], l_max_rho[b], ".-", color=color[b])
        ax3.axhline(0, color="gray", lw=0.5)
    end
    ax3.legend()

    # amp_b1 = sqrt.(ifelse.((sign(Z).*(eps_c[1] .- epsilons_b1) .> 0), sign(Z).*(eps_c[1] .- epsilons_b1), NaN)).*amp2(P.beta_a1, P.beta_a2, 2)
    # amp_b2 = sqrt.(ifelse.((sign(Z).*(eps_c[2] .- epsilons_b2) .> 0), sign(Z).*(eps_c[2] .- epsilons_b2), NaN)).*amp2(P.beta_a1, P.beta_a2, 3)

    # for s in [-1, 1]
        # ax3.plot(epsilons_b1, 1 .+ s*amp_b1, ls="dotted", color="black", lw=0.5)
        # ax3.plot(epsilons_b2, 1 .+ s*amp_b2, ls="dotted", color="black", lw=0.5)
    # end


    fig2 = PyPlot.figure(dpi=50)
    if only_branch > 0
        ax4 = PyPlot.subplot(111)
        ax4.axis("off")
    else
        ax4 = PyPlot.subplot(121)
        # ax5 = PyPlot.subplot(247)

        ax6 = PyPlot.subplot(2, 2, 2)
        ax6bis = ax6.twinx()
        ax7 = PyPlot.subplot(2, 2, 4)
        # ax8 = PyPlot.subplot(4, 4, 12)
        # ax8bis = ax8.twinx()
        # ax9 = PyPlot.subplot(4, 4, 16)
    end

    findnearest(A::Vector{Float64},t::Float64) = findmin(abs.(A.-t))[2]

    function plot_curve(idx::Tuple; init::Bool=false)
        selected_X = [l_Xs[b][i] for (b,i) in enumerate(idx)]
        for (b,X) in enumerate(selected_X)
            c = Bend.X2candidate(P, X)

            for i in 3:P.N+1
                xy[i,:] = xy[i-1,:] + P.Δs*[cos(c.θ[i-2]); sin(c.θ[i-2])]
            end

            bc = sum(xy, dims=1)/P.N
            xy .= xy .- bc

            beta = Bend.compute_beta(P, c.ρ)

            θ_dot = Bend.compute_centered_fd_θ(P, matrices, c.θ)

            (λM, λx, λy) = Bend.candidate_multipliers(P, X, matrices)

            if only_branch > 0
                color[b] = "black"
                if b != only_branch
                    continue
                end
            end
            if init
                ax4.scatter(xy[:,1], xy[:,2], s=dot_size_factor*c.ρ.^4, color=color[b])
                ax4.plot(xy[:,1], xy[:,2], lw=0.1, color=color[b])
                # ax4.plot(xy[:,1], xy[:,2], lw=1, color="white")
                ax4.set_aspect("equal")
                ax4.set_xlim(-1.5, 1.5)
                ax4.set_ylim(-1.5, 1.5)
                # ax[k].set_title("2D visualisation")

                if only_branch > 0
                    continue
                    ax4.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="off", right="off", left="off", labelleft="off")
                end

                ax6.plot(t, c.ρ, color=color[b])
                ax6.axhline(0, color="black", lw=0.5)
                ax6bis.plot(t, beta, ls="dotted", color=color[b])
                ax6.set_ylabel("ρ")
                ax6bis.set_ylabel("β(ρ)")

                ax7.plot(t, θ_dot, color=color[b])
                ax7.set_ylabel("d/ds θ")
                ax7.set_ylim(0, 10)
            else
                if only_branch > 0
                    ax4.lines[1].set_data(xy[:,1], xy[:,2])
                    ax4.collections[1].set_offsets(copy(xy))
                    ax4.collections[1].set_sizes(dot_size_factor*c.ρ)
                    ax4.tick_params(axis="both", which="both", bottom="off", top="off", labelbottom="off", right="off", left="off", labelleft="off")
                else
                    ax4.lines[b].set_data(xy[:,1], xy[:,2])
                    ax4.collections[b].set_offsets(copy(xy))
                    ax4.collections[b].set_sizes(dot_size_factor*c.ρ)
                end

                if only_branch > 0
                    continue
                end

                ax6.lines[2*b-1].set_data(t, c.ρ)
                ax6bis.lines[b].set_data(t, beta)
                ax7.lines[b].set_data(t, θ_dot)

                ax6.relim()
                ax6.autoscale_view(true,true,true)
                ax6bis.relim()
                ax6bis.autoscale_view(true,true,true)
            end
        end
        PyPlot.draw()
    end

    plot_curve(Tuple(1 for i in 1:n_branches); init=true)
    visible = 0

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
            idx = Tuple(findnearest(l_epsilons[b], event.xdata) for b in 1:maximum(keys(l_epsilons)))
            plot_curve(idx)
        end
    end

    function mouseup2(event)
        if event.inaxes == ax4
            visible = (visible + 1) % (n_branches+1)
            if visible == 0
                for (i,l) in enumerate(ax4.lines)
                    l.set_alpha(1)
                    ax4.collections[i].set_alpha(1)
                end
            else
                for (i, l) in enumerate(ax4.lines)
                    if visible == i
                        l.set_alpha(1)
                        ax4.collections[i].set_alpha(1)
                    else
                        l.set_alpha(0)
                        ax4.collections[i].set_alpha(1)
                    end
                end
            end
            PyPlot.draw()
        end
    end

    fig.canvas.mpl_connect("motion_notify_event", mousemove)
    fig.canvas.mpl_connect("button_release_event", mouseup)
    fig2.canvas.mpl_connect("button_release_event", mouseup2)

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
# PyPlot.savefig("2.pdf")

# plotall(3)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
# PyPlot.savefig("3.pdf")

# plotall(4)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
# PyPlot.savefig("4.pdf")

# plotall(5, only_branch=1)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
# PyPlot.savefig("5.pdf")

#plotall(7, only_branch=2)
#mng = PyPlot.get_current_fig_manager()
#mng.window.showMaximized()
#sleep(1)
#PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
# PyPlot.savefig("7.pdf")

# plotall(8, only_branch=2)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
# PyPlot.savefig("8.pdf")

# plotall(9)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
#
# plotall(10, only_branch=2)
# mng = PyPlot.get_current_fig_manager()
# mng.window.showMaximized()
# sleep(1)
# PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

 # plotall(100)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

 # plotall(201, dot_size_factor=10.0)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

 # plotall(204, dot_size_factor=3.0)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
 #
 # plotall(201, dot_size_factor=3.0, only_branch=3)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
 # plotall(202, dot_size_factor=3.0, only_branch=3)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
 # plotall(203, dot_size_factor=3.0, only_branch=3)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # #
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])
 # plotall(204, dot_size_factor=10.0)
 # mng = PyPlot.get_current_fig_manager()
 # mng.window.showMaximized()
 # sleep(1)
 # PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

 plotall(205, dot_size_factor=40.0)
 mng = PyPlot.get_current_fig_manager()
 mng.window.showMaximized()
 sleep(1)
 PyPlot.tight_layout(rect=[0, 0.03, 1, 0.95])

PyPlot.show()
