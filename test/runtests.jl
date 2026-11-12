using LinearAlgebra
using Random
using Flows
using CLVs
using Test

Random.seed!(0)

# ---------------------------------------------------------------------------- #
# HELPERS

# The three flows that `computeCLV` takes, for states in ℝⁿ and `N` vectors:
# the flow of the nonlinear equations `f`, and those of the tangent and adjoint
# equations `ftan` and `fadj` along a stored trajectory, all with RK4 steps of
# size `Δt`.
function flows(f, ftan, fadj, n, N, Δt)
    z   = couple(ntuple(i -> zeros(n), N)...)
    dep = CallDependency(ntuple(i -> (i,), N)...)

    flow_nl  = flow(f, RK4(zeros(n)), TimeStepConstant(Δt))
    flow_tan = flow(couple(ntuple(i -> ftan, N)...),
                    dep,
                    RK4(z, ContinuousMode(false)),
                    TimeStepFromStorage(Δt))
    flow_adj = flow(couple(ntuple(i -> fadj, N)...),
                    dep,
                    RK4(z, ContinuousMode(true)),
                    TimeStepFromStorage(Δt))

    return flow_nl, flow_tan, flow_adj
end

# `computeCLV` prints the time of each stage of the calculation. Run `f`
# without this output.
quiet(f) = redirect_stdout(f, devnull)

# Distance between the vector `a` and the unit vector `b`, once `a` is scaled
# to unit length and given the orientation of `b`.
function misalignment(a, b)
    c = dot(a, b)
    return norm(a ./ norm(a) .- sign(c) .* b)
end

# ---------------------------------------------------------------------------- #
# TESTS

include("test_angles.jl")
include("test_linear.jl")
include("test_lorenz.jl")
