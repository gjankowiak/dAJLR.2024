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

    t_min = nothing
    t_max = nothing
    e_min = nothing
    last_idx = length(t_i)

    configs = Dict([

                   "fig8" => (snapshots_idx = [25, 45],
                                 legend_position = :lc,
                                 last_idx = length(t_i),

                                 t_min = 5e-4,
                                 t_max = 5e2,
                                 e_min = 1e-6,

                                 plot_offset_x = -3,
                                 plot_scale_x = 795,
                                 plot_range = 1.75,

                                 iw = 125, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_c, \hat{\rho}_c)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

                   ])
                    "fig9" => (snapshots_idx = [3, 6, 14, 31, 44],
                                 legend_position = :rt,
                                 last_idx = length(t_i),

                                 t_min = 1e-4,
                                 t_max = 1e3,
                                 e_min = 1e-1,
                                 e_max = 2e3,

                                 rho_marker_factor = 1.5,

                                 plot_offset_x = 75,
                                 plot_scale_x = 720,
                                 plot_range = 1.5,

                                 iw = 105, #inset width

                                 yticks = ([1e-1, 1, energy_circle, 10, 100, 1000, e[1]], [M.L"10^{-1}", M.L"10^0", M.L"\hat{\mathcal{E}}_\mu(\hat{\theta}_c, \hat{\rho}_c)", M.L"10^1", M.L"10^2", M.L"10^3", M.L"\hat{\mathcal{E}}_\mu(t=0)"])),

    try
        config = configs[basename(dirname(dir * "/"))]
    catch
        println("key not found?")
        return
    end

    fig = M.Figure(size=(810, 607))

    xscale = get(config, :xscale, log10)
    yscale = get(config, :yscale, log10)

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

        xy, markersize, color = compute_xy(P, IP, s_X, matrices; color_factor=10.0, markersize_factor=markersize_factor)

        M.scatterlines!(s_ax, xy[:,1], xy[:,2], markersize=markersize, markercolor=color, linewidth=0.5, color=:gray65)

        if i > 0
            M.vlines!(ax, s_t, linewidth=1, color=:gray50)
        else
            M.arrows!(s_ax, [-1.5], [0.0], [-0.9], [0.0]; linewidth=0.5)
        end
    end

    if energy_circle > 1e-12
        M.hlines!(ax, energy_circle, linewidth=1, color=:gray20, label=M.L"\hat{\mathcal{E}}_{\mu}(\hat{\theta}_c, \hat{\rho}_c)")
    else
        M.text!(ax, M.L"\hat{\mathcal{E}}_{\mu}(\hat{\theta}_c, \mu_0) = 0 \; \rightarrow"; position=(config.t_max*0.65, config.e_min*1.5e0), align=(:right, :center), rotation=-π/2)
    end

    M.lines!(ax, t_i[2:last_idx], e[2:last_idx],                           linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}_\mu")
    M.lines!(ax, t_i[2:last_idx], energy_θ_i[2:last_idx], linestyle=:dash, linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}^{{\theta}}")
    M.lines!(ax, t_i[2:last_idx], energy_ρ_i[2:last_idx], linestyle=:dot,  linewidth=2.0, color="black", label=M.L"\hat{\mathcal{E}}_{\mu}^{{\rho}}")

    M.xlims!(ax, (config.t_min, config.t_max))
    M.ylims!(ax, (get(config, :e_min, nothing), get(config, :e_max, nothing)))

    M.axislegend(ax, position=config.legend_position)

    M.save(energy_plot_fn, fig)

    close(h5)
end
