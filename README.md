# Numerical gradient flow described in dAJLR

## Installation

- Install julia
- Clone this repository: `git clone https://git.oknaj.eu/gjankowiak/dAJLR && cd dAJLR`
- Start a new project: `julia --project=.`
- Get the packages, the will also pull and build all dependencies:
    
    ]
    add https://git.oknaj.eu/gjankowiak/EvenParam.jl
    add https://git.oknaj.eu/gjankowiak/JankoUtils.jl
    add https://git.oknaj.eu/gjankowiak/ModulatedCurves.jl

- Edit `flow_stability.jl` and run it from the julia prompt with `include("flow_stability.jl")`

## TODO

- Set options from the command line or script
