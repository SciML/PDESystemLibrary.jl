using PDESystemLibrary
PSL = PDESystemLibrary

using NeuralPDE, Lux, OptimizationOptimJL
import ModelingToolkit: Interval

for ex in PSL.all_systems
    @testset "Example: $(ex.name) with MethodOfLines.jl" begin
        dim = length(ex.ivs) # number of dimensions
        N = 8*dim
        chain = Lux.Chain(Dense(dim, N, Lux.σ), Dense(N, N, Lux.σ), Dense(N, 1))

        # Discretization
        dx = 0.05
        discretization = PhysicsInformedNN(chain, GridTraining(dx))

        @named pde_system = PDESystem(eq, bcs, domains, [x, y], [u(x, y)])
        prob = discretize(pde_system, discretization)

        #Optimizer
        opt = OptimizationOptimJL.BFGS()

        #Callback function
        callback = function (p, l)
            println("Current loss is: $l")
            return false
        end

        res = Optimization.solve(prob, opt, callback = callback, maxiters = 1000)
        phi = discretization.phi
    end
end