using ModulatedCurves
import HDF5

func_defs_fn = "function_definitions.jl"
default_func_defs_fn = string(func_defs_fn, ".default")

config_fn = "config.params.toml"
default_config_fn = string("config.params.toml", ".default")

if !isfile(func_defs_fn)
    @info "Copying '$func_defs_fn' from defaults, you can now edit it."
    cp(default_func_defs_fn, func_defs_fn)
end
include(func_defs_fn)

if !isfile(config_fn)
    @info "Copying '$config_fn' from defaults, you can now edit it." 
    cp(default_config_fn, config_fn)
end

function main(configuration_fn::String=config_fn)
    P, S, Xinit, solver_params, options = ModulatedCurves.parse_configuration(configuration_fn, function_defs)

    # HACK
    options_dict = Dict(options)

    if haskey(options_dict, :do_load)
        h5 = HDF5.h5open(options_dict[:load_snapshot_fn], "r")
        idx = options_dict[:load_snapshot_n]
        if idx > 0
            Xinit .= read(h5["s$(options_dict[:load_snapshot_n])X"])
        else
            last_idx = read(h5["s_last_idx"])
            Xinit .= read(h5["s$(s_last_idx)X"])
        end
        close(h5)
    end

    if haskey(options_dict, :output_dir)
        output_dir = options_dict[:output_dir]
        mkpath(output_dir)
        rm(joinpath(output_dir, "snapshots"); force=true, recursive=true)
        mkpath(joinpath(output_dir, "snapshots"))
        cp(configuration_fn, joinpath(output_dir, "config.toml"), force=true)
        cp(func_defs_fn, joinpath(output_dir, "functions_defs.toml"), force=true)
    end

    # check differentials
    # check_differential(P, S, Xinit)
    # return

    do_flow(P, S, Xinit, solver_params; options...)
end

# main()
