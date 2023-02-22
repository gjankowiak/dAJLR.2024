import HDF5
using UnicodePlots
import TOML
import CairoMakie
import ModulatedCurves
M = CairoMakie

include("flow.jl")

function show_monotonicity(dir::String)
    h5_fn = joinpath(dir, "data.hdf5")
    energy_plot_fn = joinpath(dir, "energy.pdf")
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

    e = energy_ρ_i + energy_θ_i

    @show P.µ

    @show energy_ρ_i[1]
    @show energy_θ_i[1]

    p = UnicodePlots.lineplot(t_i, energy_ρ_i, title="E_ρ")
    println(UnicodePlots.hline!(p, P.μ * 1.5708))

    p = UnicodePlots.lineplot(t_i, energy_θ_i, title="E_θ")
    println(UnicodePlots.hline!(p, 2.74889))

    println(UnicodePlots.lineplot(t_i, energy_θ_i + energy_ρ_i, title="E_ρ + E_θ"))

    close(h5)
end

function r(prefix; force=false)
    toml_fn = "non_monotonicity/non_monotonicity_mu=$(prefix).toml"
    data_prefix = "test/non_monotonicity_mu=$(prefix)"
    data_fn = joinpath(data_prefix, "data.hdf5")

    if force || (stat(toml_fn).mtime > stat(data_fn).mtime)
        flow(toml_fn)
    end

    show_monotonicity(data_prefix)
end
