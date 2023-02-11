# PDESystemLibrary.jl
A library of systems of partial differential equations, as defined with ModelingToolkit.jl in Julia.

This library contains various lists of systems dominated by different properties, as well as a master list.
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
    heat_1d

    1D heat equation with Dirichlet boundary conditions.
    This models the temperature of a rod over time, where the ends are held at a constant temperature.

    It is initialized with a sinusoidal profile.
    
    The equation is given by:
    ```math
    \\frac{\\partial u}{\\partial t} = D \\frac{\\partial^2 u}{\\partial x^2}
    ```
"""
heat_1d = begin
    @variables x t
    @parameters D

    Dxx = Differential(x)
    Dt = Differential(t)

    eqs = [Dt(u(t, x)) ~ D * Dxx(u(t, x))]
    bcs = [u(0, x) ~ sin(2pi * x), 
           u(t, 0) ~ 0.0, 
           u(t, 1) ~ 0.0]

    domains = [t ∈ IntervalDomain(0.0, 1.0), 
               x ∈ IntervalDomain(0.0, 1.0)]

    function u_exact(t, x; ps = (1.0,))
        exp(-4 * (pi^2) * ps[D] * t) * sin(2pi * x)
    end

    analytic = [u(t, x) => u_exact]

    @named heat_1d = PDESystem(eqs, bcs, domains, [t, x], [u(t, x)], [D => 1.0], analytic = analytic)

    heat_1d
end

# Add to the lists
push!(all_systems, heat_1d)
push!(linear_diffusion_systems, heat_1d)
```

## A note on analytic solutions
Analytic solutions should have the same argument signature as their parent variable, 
and must take a keyword argument `ps`, which may be unused.

Analytic solutions are optional, if you do not know the analytic solution do not use the keyword argument.
Downstream packages must check whether analytic solutions are present in a system before using it to handle this case.
