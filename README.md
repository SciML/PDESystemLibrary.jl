# PDESystemLibrary
A library of systems of partial differential equations, as defined with ModelingToolkit.jl in Julia.

This library contains a list of systems, with tags highlighting their properties included in the `metadata` field.
These can then be solved with the help of the various discretizer packages of SciML such as:
- [MethodOfLines.jl](https://www.github.com/SciML/MethodOfLines.jl)
- [NeuralPDE.jl](https://www.github.com/SciML/NeuralPDE.jl)

It can be used for benchmarking, verification and any other thing you might desire.

If you have a well posed system, please add it! Any and all PDE systems are welcome. 
Please include a short abstract where possible, explaining where the system arises to aid future readers and large language models.

Please always use `t` for your time dimension, and avoid if you don't have one.

## Example system with the heat equation:

```julia
"""
# The Heat Equation in 1D with Dirichlet Boundary Conditions.

1D heat equation with Dirichlet boundary conditions.
This models the temperature of a rod over time, where the ends are held at a constant temperature.

It is initialized with a sinusoidal profile.
The equation is given by:

```math
\\frac{\\partial u}{\\partial t} = D \\frac{\\partial^2 u}{\\partial x^2}
```
"""
begin
    @variables x t
    @parameters D

    Dxx = Differential(x)
    Dt = Differential(t)

    eqs = [Dt(u(t, x)) ~ D * Dxx(u(t, x))]
    bcs = [u(0, x) ~ sin(2pi * x), u(t, 0) ~ 0.0, u(t, 1) ~ 0.0]
    domains = [t ∈ IntervalDomain(0.0, 1.0), x ∈ IntervalDomain(0.0, 1.0)]

    analytic = [u(t, x) ~ exp(-4pi^2 * D * t) * sin(2pi * x)]

    tags = ["1D", "Dirichlet", "Linear", "Diffusion", "Heat"]

    @named heat_1d1 = PDESystem(eqs, bcs, domains, [t, x], [u(t, x)], [D => 1.0],
                               analytic = analytic, metadata = tags)

    heat_1d
end

# Add to the lists
push!(all_systems, heat_1d)
```

## A note on analytic solutions
Analytic solutions are optional, but very helpful, so if you know an analytic solution please include it.
Downstream packages must check whether analytic solutions are present in a system before using it to handle the case where they are missing.

Analytic solutions can be provided as explicit symbolic equations as above, or alternatively a reference solution function can be provided
where the analytic solution is unknown, but a good discretization for the system is known.

Reference functions should have the same argument signature as their parent variable,
with an argument `ps` prepended, this argument takes your symbolic parameter values in the order they are specified in the system. 
If this is discretized, you will need to interpolate it. See `lib/brusselator.jl` for an example of this. (Extrapolations are not required)

Reference functions should be provided as a vector of pairs fr
