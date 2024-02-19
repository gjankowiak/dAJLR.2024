# Numerical gradient flow described in [dAJLR][preprint]

## Installation

- Install [Julia](https://julialang.org/)
- Clone this repository: `git clone https://github.com/gjankowiak/dAJLR && cd dAJLR`
- Start a new project: `julia --project=.`
- Get the packages, the will also pull all dependencies:
```
    ] add https://git.oknaj.eu/gjankowiak/EvenParam.jl add https://git.oknaj.eu/gjankowiak/JankoUtils.jl add https://git.oknaj.eu/gjankowiak/ModulatedCurves.jl
```

- Try with the default set of parameters:
```
    include("flow.jl")
```

- To run the code, you can call the `flow` function by passing the path to a configuration file. For some of the figures, configuration files are in the `figures` directory. For examples, the solution in Figure 11 (a) is obtained with:
```
    flow("figures/fig11a.toml")
```
The results will be output to the `figures/results` directory.

- You can edit the newly created `config.params.toml` and `function_definitions.jl` run it again!

[preprint]: https://arxiv.org/abs/2308.01151
