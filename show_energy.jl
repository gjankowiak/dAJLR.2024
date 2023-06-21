import HDF5
using UnicodePlots
import TOML
import CairoMakie
import ModulatedCurves
M = CairoMakie

function compute_xy(P, IP, X, matrices; color_curvature::Bool=false, color_factor::Float64=10.0, markersize_factor::Real=1.0)
    c = ModulatedCurves.X2candidate(P, X)
    xy = zeros(P.N+1,2)
    for i in 2:P.N+1
        xy[i,:] = xy[i-1,:] + IP.Δs*[cos(c.θ[i-1]); sin(c.θ[i-1])]
    end

    bc = sum(xy, dims=1)/P.N
    xy = xy .- bc

    θ_dot = ModulatedCurves.compute_centered_fd_θ(P, matrices, c.θ)

    if !color_curvature
        # mapped_val = r -> clamp(abs(r)/color_factor, 0, 1.0)
        # colors = map(r -> (r > 0 ? M.Colors.RGBA(0.0, 0.0, 255.0, mapped_val(r)) : M.Colors.RGBA(255.0, 0.0, 0.0, mapped_val(r))), [c.ρ; c.ρ[1]])
        # markersize = 3*markersize_factor
        zipped_ρ_θ_dot = [[(ρ, θ_dot[i]) for (i,ρ) in enumerate(c.ρ)]; (c.ρ[1], θ_dot[1])]
        colors = map(t -> t[1] > 0 ? "blue" : "red", zipped_ρ_θ_dot)

        markersize = markersize_factor*abs.([c.ρ; c.ρ[1]])
    else
        zipped_ρ_θ_dot = [[(ρ, θ_dot[i]) for (i,ρ) in enumerate(c.ρ)]; (c.ρ[1], θ_dot[1])]
        colors = map(t -> t[1] > 0 ? (t[2] >= 0.0 ? "blue" : "lightblue") : (t[2] >= 0.0 ? "red" : "pink"), zipped_ρ_θ_dot)

        markersize = markersize_factor*abs.([c.ρ; c.ρ[1]])
    end

    return xy, markersize, colors
end

function format_ticks(values)
    map(values) do v
        e = round(log10(v); digits=1)
        try
            e = Int(e)
        catch
        end
        return M.L"10^{%$(e)}"
    end
end

function show_energy(dir::String)
    h5_fn = joinpath(dir, "data.hdf5")
    bn = basename(dir)
    if bn == ""
        bn = basename(dirname(dir))
    end
    @info dir, bn
    energy_plot_fn = joinpath(dir, "energy_$(basename(dir)).pdf")
    toml_fn = joinpath(dir, "config.toml")

    config = TOML.parsefile(toml_fn)

    P = ModulatedCurves.params_from_toml(config)
    IP = ModulatedCurves._compute_intermediate(P)
    @show IP.Δs


    h5 = HDF5.h5open(h5_fn, "r")

    n_last = read(h5["n_last"])

    energy_ρ_i = read(h5["energy_rho"])[1:n_last-1]
    energy_θ_i = read(h5["energy_theta"])[1:n_last-1]
    t_i = read(h5["t"])[1:n_last-1]

    energy_circle = read(h5["energy_circle"])

    Xinit = read(h5["s0X"])
    w = ModulatedCurves.compute_winding_number(P, Xinit)
    matrices = ModulatedCurves.assemble_fd_matrices(P, IP; winding_number=w)

    e = energy_ρ_i + energy_θ_i
    # plt = lineplot(t_i .+ 1, e, width = :auto, height = :auto, xscale=:log10)
    # hline!(plt, energy_circle)
    # 
    # println(plt)

    # plt = lineplot(t_i .+ 1, energy_ρ_i, width = :auto, height = :auto, xscale=:log10)
    # hline!(plt, 0.0)
    # println(plt)

    t_min = nothing
    t_max = nothing
    e_min = nothing
    last_idx = length(t_i)

    configs = Dict(
                   # ["mu=3e-3_rhoA=2.0_rA=0.4" => (snapshots_idx = [3, 6, 14, 22, 36, 44, 53], # 11, 40
                   ["mu=3e-3_rhoA=2.0_rA=0.4" => (snapshots_idx = [3, 6, 14, 31, 44],
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 1e3,
                                 e_min = 1e-1,

                                 plot_offset_x = 90,
                                 plot_scale_x = 690,
                                 plot_range = 1.5,

                                 iw = 75, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                    "mu=1e-3_rhoA=2.0_rA=0.4" => (snapshots_idx = [3, 6, 14, 31, 44],
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 1e3,
                                 e_min = 1e-1,

                                 plot_offset_x = 90,
                                 plot_scale_x = 690,
                                 plot_range = 1.5,

                                 iw = 75, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                    "mu=1e-1_rhoA=2.0_rA=0.4" => (snapshots_idx = [3, 7, 12, 15, 23, 32, 38, 54],
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 1e3,
                                 e_min = 1e-3,

                                 plot_offset_x = 90,
                                 plot_scale_x = 695,
                                 plot_range = 1.5,

                                 iw = 45, #inset width

                                 yticks = ([1e-2, 1, energy_circle, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                    "mu=1e1_rhoA=2.0_rA=0.4" => (snapshots_idx = [3, 7, 12, 15, 23],
                                 legend_position = :lb,
                                 last_idx = findfirst(x -> x > 5e-1, t_i),
                                 t_min = 1e-4,
                                 t_max = 5e-1,
                                 e_min = 1e-5,

                                 iw = 80, #inset width,

                                 plot_offset_x = 90,
                                 plot_scale_x = 1125,
                                 plot_range = 1.2,

                                 yticks = ([1e-2, 1, energy_circle, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                   # "mu=1e-3_rhoA=0.01_rA=0.01" => (snapshots_idx = [22, 37, 48, 54], # 40
                   "mu=1e-3_rhoA=0.01_rA=0.01" => (snapshots_idx = [22, 48],
                                 legend_position = :lc,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 4e3,
                                 e_min = 1e-6,

                                 plot_offset_x = 90,
                                 plot_scale_x = 640,
                                 plot_range = 1.5,

                                 iw = 85, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                   "mu=1e-3_rhoA=0_rA=0.01" => (snapshots_idx = [22, 38, 56], # 40
                                 legend_position = :lc,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 4e3,
                                 e_min = 1e-4,

                                 plot_offset_x = 90,
                                 plot_scale_x = 725,
                                 plot_range = 1.75,

                                 iw = 125, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                   "mu=3e-3_rhoA=2.0_rA=0.4_c0=1" => (snapshots_idx = [4, 14, 21, 50], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 4e3,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 plot_offset_x = 90,
                                 plot_scale_x = 638,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 85, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, energy_circle, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),

                                 title = M.L"c_0 = 1,\; \mu = 3 \times 10^{-3}"
                                ),

                   "c0=3.mu_low.roundinit" => (snapshots_idx = [4, 7, 10, 17], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-2,
                                 t_max = 1e3,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 830,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 110, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, energy_circle, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_0, \hat{\rho}_0)", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "mu=1e-2_rhoA=0.01_rA=0.01.omega=2" => (snapshots_idx = [4, 14, 20, 21, 43], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-2,
                                 t_max = 1e4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=0.1" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 1e-3,
                                 e_min = 1e-3,
                                 e_max = 1e-2,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=1.5" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=high" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=low" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=low2" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "non_monotonicity_mu=med" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   "a=5e-1" => (snapshots_idx = [1, 2, 3], # 40
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-3,
                                 t_max = 1e-4,
                                 e_min = 1e-5,
                                 e_max = 1e4,

                                 rho_marker_factor = 0.5,

                                 plot_offset_x = 90,
                                 plot_scale_x = 700,
                                 plot_range = 1.75,

                                 xscale = log10,
                                 yscale = log10,

                                 iw = 70, #inset width

                                 # yticks = M.WilkinsonTicks(5))
                                 yticks = ([1e-2, 1, 100, e[1]], [M.L"10^{-2}", M.L"10^0", M.L"10^2", M.L"\hat{\mathcal{E}}_\mu(t=0)"]),
                                ),
                   ])

    try
        config = configs[basename(dirname(dir * "/"))]
    catch
        println("key not found?")
        return
    end

    fig = M.Figure()

    xscale = get(config, :xscale, log10)
    yscale = get(config, :yscale, log10)

    @show P.µ

    # xscale = M.Makie.pseudolog10
    ax = M.Axis(fig[2:5,1:10], xscale=xscale, yscale=yscale, xlabel=M.L"$t$", yticks=config.yticks, xtickformat=format_ticks, title=get(config, :title, M.L"\mu = %$(P.µ)"))

    s_t_i = []
    for (n,i) in enumerate(config.snapshots_idx)
        try
            push!(s_t_i, log10(read(h5["s$(i)t"])))
        catch
            println(backtrace())
            println("Error")
            close(h5)
            return
        end
    end

    if !haskey(config, :t_min)
        log_t_min_rounded  = round(log10(t_i[2]), RoundDown)
    else
        log_t_min_rounded  = round(log10(config.t_min), RoundDown)
    end
    log_t_max_rounded  = round(log10(t_i[end]), RoundUp)
    log_t_delta = log_t_max_rounded - log_t_min_rounded

    left = log_t::Float64 -> (log_t-log_t_min_rounded)/log_t_delta*config.plot_scale_x + config.plot_offset_x

    bboxes = [M.BBox(left(log_t)-config.iw/2, left(log_t) + config.iw/2, 490, 490+1config.iw) for log_t in s_t_i]
    limits = [(-config.plot_range, config.plot_range, -config.plot_range, config.plot_range) for log_t in s_t_i]

    pushfirst!(config.snapshots_idx, 0)
    pushfirst!(bboxes, M.BBox(10, config.iw+40, 490, 490+1config.iw))
    pushfirst!(limits, (-config.plot_range-1.5, config.plot_range, -config.plot_range, config.plot_range))

    snapshots_axes = [M.Axis(fig, bbox=bboxes[n],  autolimitaspect=1.0, limits=limits[n]) for (n,i) in enumerate(config.snapshots_idx)]
    # snapshots_axes = [M.Axis(fig, bbox=bboxes[n],  aspect=M.AxisAspect(1), limits=limits[n], backgroundcolor=:blue) for (n,i) in enumerate(snapshots_idx)]

    # offset = 10t_i[2]
    # M.scatterlines!(ax, offset .+ t_i, e, color="black", label=M.L"\hat{\mathcal{E}}_\mu")
    # M.lines!(ax, offset .+ t_i, e, color="black", label=M.L"\hat{\mathcal{E}}_\mu")
    # M.lines!(ax, offset .+ t_i, linestyle=:dash, color="black", energy_θ_i, label=M.L"\hat{\mathcal{E}}_1")
    # M.lines!(ax, offset .+ t_i, linestyle=:dot, color="black", energy_ρ_i, label=M.L"\hat{\mathcal{E}}_2")

    if haskey(config, :annotations)
        M.annotations!(ax, config.annotations.text, config.annotations.points)
    end

    markersize_factor = get(config, :rho_marker_factor, 1.0)

    for (n,i) in enumerate(config.snapshots_idx)
        s_X = read(h5["s$(i)X"])
        s_t = read(h5["s$(i)t"])
        s_n = read(h5["s$(i)n"])
        s_ax = snapshots_axes[n]
        M.hidedecorations!(s_ax)
        M.hidespines!(s_ax)

        #@show s_n, s_t

        xy, markersize, color = compute_xy(P, IP, s_X, matrices; color_factor=10.0, markersize_factor=markersize_factor)

        M.scatterlines!(s_ax, xy[:,1], xy[:,2], markersize=markersize, markercolor=color, linewidth=0.5, color=:gray65)

        if i > 0
            M.vlines!(ax, s_t, linewidth=1, color=:gray50)
        else
            M.arrows!(s_ax, [-1.5], [0.0], [-0.9], [0.0]; linewidth=0.5)
        end
    end

    M.hlines!(ax, e[1], linewidth=1, color=:gray20)
    if energy_circle > 1e-12
        M.hlines!(ax, energy_circle, linewidth=1, color=:gray20, label=M.L"\hat{\mathcal{E}}_{\mu}(\hat{\theta}_0, \hat{\rho}_0)")
    else
        # M.arrows!(ax, [config.t_max*0.7], [config.e_min*1.4e1], [0.0], [-1.3e-4]; linewidth=0.5)
        M.text!(ax, M.L"\hat{\mathcal{E}}_{\mu}(\hat{\theta}_0, \mu_0) = 0 \; \rightarrow"; position=(config.t_max*0.65, config.e_min*1.5e0), align=(:right, :center), rotation=-π/2)
    end

    @show energy_ρ_i[1]
    @show energy_θ_i[1]

    p = UnicodePlots.lineplot(t_i, energy_ρ_i, title="E_ρ")
    println(UnicodePlots.hline!(p, P.μ * 1.5708))

    p = UnicodePlots.lineplot(t_i, energy_θ_i, title="E_θ")
    println(UnicodePlots.hline!(p, 2.74889))

    println(UnicodePlots.lineplot(t_i, energy_θ_i + energy_ρ_i, title="E_ρ + E_θ"))

    M.lines!(ax, t_i[2:last_idx], e[2:last_idx],                           linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}_\mu")
    M.lines!(ax, t_i[2:last_idx], energy_θ_i[2:last_idx], linestyle=:dash, linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}^{{\theta}}")
    M.lines!(ax, t_i[2:last_idx], energy_ρ_i[2:last_idx], linestyle=:dot,  linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}_{\mu}^{{\rho}}")

    M.xlims!(ax, (config.t_min, config.t_max))
    M.ylims!(ax, (get(config, :e_min, nothing), get(config, :e_max, nothing)))

    M.axislegend(ax, position=config.legend_position)

    M.save(energy_plot_fn, fig)

    close(h5)
end
