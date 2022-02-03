using ModulatedCurves

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

function main()
    P, S, Xinit, solver_params, options = ModulatedCurves.parse_configuration(config_fn, function_defs)

    do_flow(P, S, Xinit, solver_params; options...)
end

main()
