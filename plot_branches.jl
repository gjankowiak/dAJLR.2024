push!(LOAD_PATH, "src")

import Bend
import PyPlot
import Serialization
import Interpolations

# branch 1
epsilons_b1 = [range(0.35, 0.21, length=50); 2e-1; 1e-1; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2; 5*10. .^range(-2, -4; length=10)[2:6]]
energy_b1 = zeros(size(epsilons_b1))
amp_b1 = zeros(size(epsilons_b1))

eps_c = [0.353553390593273; 0.23570226039551584]

# branch 2
epsilons_b2 = [0.2355; 0.23545; 0.23525; range(0.235, 0.2, step=-0.0005); range(0.2, 0.1, length=10)[2:end]; 9e-2; 8e-2; 7e-2; 6e-2; 5e-2; 5*10. .^range(-2, -4; length=10)[2:6]]
energy_b2 = zeros(size(epsilons_b2))
amp_b2 = zeros(size(epsilons_b2))

# Xs_b1 = Serialization.deserialize("results_case_1/branch_1.dat")
Xs_b1 = Serialization.deserialize("results_case_1/branch_1.dat")
P_b1 = Serialization.deserialize("results_case_1/branch_1_P.dat")
Ps_b1 = [Bend.copy(P_b1) for i in 1:length(Xs_b1)]

matrices = Bend.assemble_fd_matrices(P_b1)

P = P_b1
N = P.N
Δs = P.Δs

t = collect(range(0, 2π, length=N+1))[1:N]

xy = zeros(N+1, 2)
xy2 = zeros(N+1, 2)
xy[2,:] = [Δs 0]

for i in 1:length(Xs_b1)
    Ps_b1[i].epsilon = epsilons_b1[i]
    P = Ps_b1[i]
    energy_b1[i] = Bend.compute_energy(P, matrices, Xs_b1[i])

    c = Bend.X2candidate(P, Xs_b1[i])
    min_rho, max_rho = extrema(c.ρ)
    amp_b1[i] = max_rho - min_rho
end


# Xs_b2 = Serialization.deserialize("results_case_1/branch_2.dat")
Xs_b2 = Serialization.deserialize("results_case_1/branch_2.dat")
P_b2 = Serialization.deserialize("results_case_1/branch_2_P.dat")
Ps_b2 = [Bend.copy(P_b2) for i in 1:length(Xs_b2)]

for i in 1:length(Xs_b2)
    Ps_b2[i].epsilon = epsilons_b2[i]
    P = Ps_b2[i]
    energy_b2[i] = Bend.compute_energy(P, matrices, Xs_b2[i])

    c = Bend.X2candidate(P, Xs_b2[i])
    min_rho, max_rho = extrema(c.ρ)
    amp_b2[i] = max_rho - min_rho
end

itp_energy_b1 = Interpolations.LinearInterpolation(reverse(epsilons_b1), reverse(energy_b1))

energy_diff = itp_energy_b1(epsilons_b2) - energy_b2

fig = PyPlot.figure(dpi=50)
ax1 = PyPlot.subplot(231)
ax1.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
ax1.plot(epsilons_b1, energy_b1, label="branch 1", ".-")
ax1.plot(epsilons_b2, energy_b2, label="branch 2", ".-")
ax1.axhline(π, label="circle", color="black")
for eps in eps_c
    ax1.axvline(eps, color="black", lw=0.5)
end
ax1.set_title("Energy")
ax1.legend()

ax2 = PyPlot.subplot(234, sharex=ax1)
ax2.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
ax2.plot(epsilons_b2, energy_diff, ".-")
ax2.axhline(0, color="black", lw=0.5)
ax2.axvline(0, color="black", lw=0.5)
ax2.set_title("Energy difference between branches 1 and 2")

ax3 = PyPlot.subplot(132)
ax3.axvline(epsilons_b1[1], ls="dotted", color="black", lw=0.5)
ax3.plot(epsilons_b1, amp_b1, label="branch 1", ".-")
ax3.plot(epsilons_b2, amp_b2, label="branch 2", ".-")
ax3.axhline(0, label="circle", color="black")
for eps in eps_c
    ax3.axvline(eps, color="black", lw=0.5)
end
ax3.set_title("Amplitude in rho")
ax3.legend()

ax4 = PyPlot.subplot(233)
ax5 = PyPlot.subplot(236)

findnearest(A::Vector{Float64},t::Float64) = findmin(abs.(A.-t))[2]

function plot_curve(i1::Int64, i2::Int64; init::Bool=false)
    X1 = Xs_b1[i1]
    X2 = Xs_b2[i2]
    ax = [ax4, ax5]
    color = ["#1f77b4", "#ff7f0e"]
    xy_pos = [xy, xy2]
    for (k,X) in enumerate((X1, X2))
        xy = xy_pos[k]
        c = Bend.X2candidate(P, X)

        for i in 3:N+1
            xy[i,:] = xy[i-1,:] + Δs*[cos(c.θ[i-2]); sin(c.θ[i-2])]
        end

        bc = sum(xy, dims=1)/N
        xy .= xy .- bc

        if init
            ax[k].plot(xy[:,1], xy[:,2], lw=0.1, color=color[k])
            ax[k].scatter(xy[:,1], xy[:,2], s=3*c.ρ, color=color[k])
            ax[k].set_aspect("equal")
            ax[k].set_xlim(-1.5, 1.5)
            ax[k].set_ylim(-1.5, 1.5)
            ax[k].set_title("2D visualisation")
        else
            ax[k].lines[1].set_data(xy[:,1], xy[:,2])
            ax[k].collections[1].set_offsets(xy)
            ax[k].collections[1].set_sizes(3*c.ρ)
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


PyPlot.show()
