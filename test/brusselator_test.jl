using PDESystemLibrary
using Test

# Reference values of the interpolated TRBDF2 solution, recorded from master before
# init_brusselator_2d stopped boxing N in brusselator_2d_loop (a runtime fix, not a numerical one).
@testset "brusselator_2d reference solution" begin
    system = PDESystemLibrary.brusselator_2d()
    analytic = Dict(string(k) => f for (k, f) in system.analytic_func)
    points = [
        (0.25, 0.5, 1.0), (0.5, 0.5, 2.5), (0.3, 0.6, 5.0),
        (0.7, 0.2, 8.0), (0.9, 0.9, 11.5), (0.123, 0.456, 3.21),
    ]
    v_expected = [
        2.8811755674441066, 4.108333882232245, 3.018906839197202,
        3.433758319348464, 0.729693635997792, 4.573621047551331,
    ]
    u_expected = [
        0.5857110711273517, 0.40724025900033406, 2.635224785941256,
        0.42123036060154495, 0.5945651848870898, 0.4518910002856509,
    ]
    for (i, point) in enumerate(points)
        @test analytic["u(x, y, t)"](nothing, point...) ≈ u_expected[i] rtol = 1.0e-8
        @test analytic["v(x, y, t)"](nothing, point...) ≈ v_expected[i] rtol = 1.0e-8
    end
end
