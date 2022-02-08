# Numerical gradient flow described in [dAJLR][preprint]

## Installation

- Install [Julia](https://julialang.org/)
- Clone this repository: `git clone https://git.oknaj.eu/gjankowiak/dAJLR && cd dAJLR`
- Start a new project: `julia --project=.`
- Get the packages, the will also pull and build all dependencies:
```
    ]
    add https://git.oknaj.eu/gjankowiak/EvenParam.jl
    add https://git.oknaj.eu/gjankowiak/JankoUtils.jl
    add https://git.oknaj.eu/gjankowiak/ModulatedCurves.jl
```

- Try with the default set of parameters:
```
    include("flow.jl")
```
- Edit the newly created `config.params.toml` and `function_definitions.jl` run it again!

[preprint]: https://www.overleaf.com/project/6082890f9f7f0bb5a15842b6
