# The Lorenz equations, with σ = 10, ρ = 28 and β = 8/3, on a chaotic
# trajectory. The CLVs and the adjoint CLVs are not known in closed form, and
# the tests check the properties that define them, with integrations that
# `computeCLV` does not use.

function lorenz(t, x, dxdt)
    dxdt[1] = 10*(x[2] - x[1])
    dxdt[2] = 28*x[1] - x[2] - x[1]*x[3]
    dxdt[3] = -8/3*x[3] + x[1]*x[2]
    return dxdt
end

# tangent equations, the Jacobian at u applied to y
function lorenz_tan(t, u, y, dydt)
    dydt[1] = 10*(y[2] - y[1])
    dydt[2] = (28 - u[3])*y[1] - y[2] - u[1]*y[3]
    dydt[3] = u[2]*y[1] + u[1]*y[2] - 8/3*y[3]
    return dydt
end

# adjoint equations, the transpose of the Jacobian at u applied to w
function lorenz_adj(t, u, w, dwdt)
    dwdt[1] = -10*w[1] + (28 - u[3])*w[2] + u[2]*w[3]
    dwdt[2] = 10*w[1] - w[2] + u[1]*w[3]
    dwdt[3] = -u[1]*w[2] - 8/3*w[3]
    return dwdt
end

@testset "Lorenz equations                       " begin

    # ---- time step, spin-up, window and sampling ----
    Δt       = 0.005
    tspin    = 25.0
    T        = 2.0
    Nsplits  = 5
    oneevery = 20

    # ---- a point on the attractor ----
    flow_nl, flow_tan, flow_adj = flows(lorenz, lorenz_tan, lorenz_adj, 3, 3, Δt)
    x0 = flow_nl([1.0, 1.0, 20.0], (0, 20))

    # ---- the three CLVs and their adjoints ----
    ts, Γs, Θs = quiet() do
        computeCLV(flow_nl, flow_tan, flow_adj, copy(x0), () -> randn(3),
                   tspin, T, Nsplits, 3, oneevery)
    end

    # ---- the trajectory over the window ----
    # The spin-up goes in the same chunks as in `computeCLV`, so that the
    # trajectory is the same to the last bit.
    x = copy(x0)
    t = 0
    for i in 1:Nsplits
        flow_nl(x, (t, t + tspin/Nsplits))
        t += tspin/Nsplits
    end
    store = RAMStorage(x)
    flow_nl(x, (tspin, tspin + T), store)

    # states and times of the trajectory at the samples of the CLVs
    xk = samples(store)[1:oneevery:end]
    tk = times(store)[1:oneevery:end]

    @test tk .- tspin ≈ ts
    @test maximum(abs(norm(Γ[j]) - 1) for Γ in Γs, j in 1:3) < 1e-12

    # ---- the two families are dual ----
    @test maximum(abs(dot(Θ[j], Γ[k]) - (j == k))
                  for (Γ, Θ) in zip(Γs, Θs), j in 1:3, k in 1:3) < 1e-10

    # ---- the neutral CLV is the direction of the flow ----
    # This holds for the exact tangent equations, hence to the error of the
    # time discretisation here, 3e-6 with this time step.
    dxdt = zeros(3)
    err  = 0.0
    for i in eachindex(ts)
        lorenz(tk[i], xk[i], dxdt)
        err = max(err, misalignment(Γs[i][2], dxdt ./ norm(dxdt)))
    end
    @test err < 1e-4

    # ---- covariance ----
    # Over a sampling interval the tangent flow maps the CLV γ onto μγ at the
    # end of the interval, and the adjoint flow maps the adjoint CLV θ at the
    # end onto μθ at the start, with the same factor μ.
    errΓ = 0.0
    errΘ = 0.0
    for i in 1:length(ts)-1
        v = copy(Γs[i])
        w = copy(Θs[i+1])
        flow_tan(v, store, (tk[i], tk[i+1]))
        flow_adj(w, store, (tk[i+1], tk[i]))
        for j in 1:3
            μ    = dot(Θs[i+1][j], v[j])
            errΓ = max(errΓ, norm(v[j] .- μ .* Γs[i+1][j]) / norm(v[j]))
            errΘ = max(errΘ, norm(w[j] .- μ .* Θs[i][j]) / norm(w[j]))
        end
    end
    @test errΓ < 1e-10
    @test errΘ < 1e-10

    # ---- truncation: the first two vectors do not depend on N ----
    flow_nl, flow_tan, flow_adj = flows(lorenz, lorenz_tan, lorenz_adj, 3, 2, Δt)
    _, Γs2, Θs2 = quiet() do
        computeCLV(flow_nl, flow_tan, flow_adj, copy(x0), () -> randn(3),
                   tspin, T, Nsplits, 2, oneevery)
    end
    errΓ = 0.0
    errΘ = 0.0
    for i in eachindex(ts), j in 1:2
        s    = sign(dot(Γs2[i][j], Γs[i][j]))
        errΓ = max(errΓ, norm(Γs2[i][j] .- s .* Γs[i][j]))
        errΘ = max(errΘ, norm(Θs2[i][j] .- s .* Θs[i][j]) / norm(Θs[i][j]))
    end
    @test errΓ < 1e-6
    @test errΘ < 1e-6
end
