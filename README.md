# Numerical gradient flow evolving heterogeneous elastic wires

As described in:

[dAJLR] Dall'Acqua, Jankowiak, Langer, Rupp, *Conservation, convergence, and computation for evolving heterogeneous elastic wires*.

See the [preprint][preprint] on the arXiv or the [published version][published].

The core code is implemented in a separate package: [ModulatedCurves.jl][mc].

## Installation

- Install [Julia](https://julialang.org/)
- Clone this repository: `git clone https://github.com/gjankowiak/dAJLR.2024 && cd dAJLR.2024`
- Start a new project: `julia --project=.`
- Get the packages, the will also pull all dependencies:
```
    ] add https://github.com/gjankowiak/ModulatedCurves.jl
```

## Usage

- Try with the default set of parameters:
```
    include("flow.jl")
```

- To run the code, you can call the `flow` function by passing the path to a configuration file. For some of the figures, configuration files are in the `figures` directory. For examples, the solution in Figure 11 (a) is obtained with:
```
    flow("figures/fig11/a.toml")
```
The results will be output to the `figures/results` directory.

- You can edit the newly created `config.params.toml` and `function_definitions.jl` run it again!

The energy plots (Figures 8 and 9) can be generated with
```
    include("show_energy")
    show_energy("figures/results/fig8")
    show_energy("figures/results/fig9")
```
The resulting `.pdf` file will be located in the corresponding result directory.


## Capture of an hand-drawn initial curve

The `capture` directory contains two scripts that can be used to capture mouse/pen table input to use as initial data for the curve. They require `opencv`, `numpy` and `matplotlib`. 
First run `python capture.py`, draw a curve with the mouse and hit `ESC` when done. Then run `python postproc.py`, which will output a `.csv` file usable as initial data
(see the `filename` key in the configuration files, i.e. `figures/fig12.toml`).

Note that in this case, the initial density $\rho$ will be constant along the curve.

[preprint]: https://arxiv.org/abs/2308.01151
[published]: https://doi.org/10.1137/23M159086X
[mc]: https://github.com/gjankowiak/ModulatedCurves.jl
