# Julia_DE_templates
This repository is a template showing how to use Julia to integrate differential equations.

template.jl shows how to solve a simple ODE problem using the DifferentialEquations.jl package.
SH_template.jl combines it with SphericalHarmonicExpansions.jl for typical ion/electron trap case.

## Environment setup
To run the code, you need to install Julia (https://julialang.org/install/) and set up the environment. You can set up the environment by opening Julia's REPL in the current folder (the folder with `Project.toml`) and running the following commands:

```julia
using Pkg
Pkg.activate(".") # activate the environment in the current directory
Pkg.instantiate() # install the packages in the environment
```
This will create a `Manifest.toml` file in the current directory, which records the exact versions of the packages used in this environment. You can then run the code in `template.jl` or `SH_template.jl` by opening them in a Julia editor (e.g., VSCode with Julia extension) and running the code.
## Usage
You can modify the parameters in the code to fit your own problem. For example, in `SH_template.jl`, you can replace the fitted polynomials `RF_p` and `DC_p` with your own fitted polynomials. You can also change the initial conditions, time span, and numerical parameters to suit your needs. After modifying the code, you can run it to see the results.
