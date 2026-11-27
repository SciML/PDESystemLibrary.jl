using PDESystemLibrary
using Test

# Samples of the interpolated TRBDF2 reference solution.
# The tolerance covers platform differences in the adaptive stiff solver (about 2e-8 relative on 32-bit).
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
        @test analytic["u(x, y, t)"](nothing, point...) ≈ u_expected[i] rtol = 1.0e-6
        @test analytic["v(x, y, t)"](nothing, point...) ≈ v_expected[i] rtol = 1.0e-6
    end
end

# Calling through a function barrier keeps the callback specialized, so any allocation measured
# is the callback's own. A captured variable boxed as Core.Box allocated ~1 MB per call
# (917 KB on 32-bit); the budget leaves room for incidental allocations on other Julia versions.
function rhs_allocations(f, du, u, p, t)
    f(du, u, p, t)
    return @allocated f(du, u, p, t)
end

@testset "brusselator_2d reference right-hand side" begin
    prob, xyd = PDESystemLibrary.brusselator_2d_reference_problem()
    f = prob.f.f
    u = fill(1.0, size(prob.u0))
    du_before = similar(u)
    du_after = similar(u)

    @test rhs_allocations(f, du_before, u, prob.p, 1.0) < 1024

    # A uniform state has no diffusion, so the forcing onset at t = 1.1 is the only difference.
    f(du_before, u, prob.p, 1.0)
    f(du_after, u, prob.p, 1.2)
    inside = (argmin(abs.(xyd .- 0.3)), argmin(abs.(xyd .- 0.6)))
    outside = (length(xyd), length(xyd))
    @test du_after[inside..., 1] - du_before[inside..., 1] ≈ 5.0
    @test du_after[outside..., 1] == du_before[outside..., 1]
    @test du_after[:, :, 2] == du_before[:, :, 2]
end
