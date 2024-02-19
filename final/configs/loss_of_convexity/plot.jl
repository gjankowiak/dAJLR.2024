import Pkg
Pkg.activate("../../..")

using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)

import TOML
import HDF5

using ModulatedCurves

using ArgParse

function_defs = nothing

function expressionize_dict(src_fn)
    dst_fn = joinpath(dirname(src_fn), replace(basename(src_fn), r".jl$" => ".expr.jl"))

    src = open(src_fn, "r")
    dst = open(dst_fn, "w")

    for l in readlines(src; keep=true)
        if contains(l, "\" =>")
            write(dst, replace(l, "=> (" => "=> :((",
                                  r",$" => "),",
                                  r"\]$" => "])",
                                  r"\)$" => "))"
                                 ))
        else
            write(dst, l)
        end
    end

    close(src)
    close(dst)

    return dst_fn
end

function generate_from_dict(d; T::DataType=Symbol)
    if T == Symbol
        return Dict{T,Dict{String,Function}}( (k, generate_from_dict(v, T=String)) for (k,v) in d )
    else
        return Dict{T,Function}( (k, @RuntimeGeneratedFunction(v)) for (k,v) in d )
    end
end

function load_functions_dict(src_fn)
    dst_fn = expressionize_dict(src_fn)
    fd = include(dst_fn)
    return generate_from_dict(fd)
end

function parse_commandline()

    s = ArgParseSettings()
    @add_arg_table! s begin
        # "--opt1"
        #     help = "an option with an argument"
        "--box-hlines"
            help = "Trace horizontal for the min/max values of y for the initial condition"
            action = :store_true
        "--circle"
            help = "Draw the unit circle"
            action = :store_true
        "--range"
            help = "Plot range"
            arg_type = Float64
            default = 1.0
        "--monochrome-density"
            help = "Use only one color for the density"
            action = :store_true
        "--rho-factor", "-r"
            help = "Factor controlling the thickness of the lines"
            arg_type = Float64
            default = 100.0
    end

    return parsed_args = parse_args(ARGS, s)
end

function plot(dir_name::String, plot_options::Dict)
    configuration_fn = joinpath(dir_name, "config.toml")
    function_defs_fn = joinpath(dir_name, "function_definitions.jl")

    generated_function_defs = load_functions_dict(function_defs_fn)

    P, S, Xinit, solver_params, options = ModulatedCurves.parse_configuration(configuration_fn, generated_function_defs)

    h5 = HDF5.h5open(joinpath(dir_name, "data.hdf5"))

    snapshots_dir = joinpath(dir_name, "snapshots")

    IP = compute_intermediate(P, S)
    fig, update_plot = init_plot(P, IP, S, Xinit; plain=true, rho_factor=plot_options["rho-factor"], plot_range=plot_options["range"], no_circle=!plot_options["circle"], monochrome=plot_options["monochrome-density"], show_nodes=true, box_hlines=plot_options["box-hlines"])

    n = read(h5["s_last_idx"])
    for i in 0:n
        X = read(h5["s$(i)X"])
        t = read(h5["s$(i)t"])
        n = read(h5["s$(i)n"])
        update_plot(X, "", nothing)

        rounded_t = round(t; digits=5)

        if i == 0
            save_path = joinpath(snapshots_dir, "n=00000_t=0.pdf")
        elseif i == n
            save_path = joinpath(snapshots_dir, "n=$(lpad(n, 5, "0"))_t=$rounded_t.last.pdf")
        else
            save_path = joinpath(snapshots_dir, "n=$(lpad(n, 5, "0"))_t=$rounded_t.pdf")
        end
        M.save(save_path, fig)
    end
    close(h5)
end

function main()
    parsed_args = parse_commandline()
    @show parsed_args
    options_dict = Dict(parsed_args...)
    plot("../../results/loss_of_convexity_1440", options_dict)
end

main()