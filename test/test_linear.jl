# The fixed point x = 0 of the linear system x' = A x. The tangent equations are
# v' = A v at every time, so that the CLVs are the right eigenvectors of A, by
# decreasing eigenvalue, at every time of the window.

@testset "fixed point of a linear system         " begin

    # ---- the matrix, from its eigenvalues and its unit eigenvectors ----
    n = 5
    λ = [0.8, 0.3, -0.1, -0.6, -1.5]
    V = Matrix(1.0I, n, n) .+ 0.4 .* [sin(i + 2j) for i in 1:n, j in 1:n]
    V = V ./ sqrt.(sum(abs2, V; dims=1))
    A = V * Diagonal(λ) / V

    # ---- nonlinear, tangent and adjoint equations ----
    f(t, x, dxdt)       = mul!(dxdt, A, x)
    ftan(t, u, v, dvdt) = mul!(dvdt, A, v)
    fadj(t, u, w, dwdt) = mul!(dwdt, A', w)

    # ---- all the vectors, and the first three only ----
    for N in (5, 3)
        flow_nl, flow_tan, flow_adj = flows(f, ftan, fadj, n, N, 0.01)

        # a window of 20 sampling intervals, after a spin-up of 60 time units
        ts, Γs = quiet() do
            computeCLV(flow_nl, flow_tan, flow_adj, zeros(n), () -> randn(n),
                       60.0, 2.0, 4, N, 10)
        end

        @test ts ≈ 0:0.1:2
        @test length(Γs) == 21

        # the CLVs are the eigenvectors, up to their orientation
        @test maximum(misalignment(Γ[j], V[:, j]) for Γ in Γs, j in 1:N) < 1e-8
        @test maximum(abs(norm(Γ[j]) - 1) for Γ in Γs, j in 1:N) < 1e-12
    end

    # ---- a window that is not a whole number of sampling intervals ----
    flow_nl, flow_tan, flow_adj = flows(f, ftan, fadj, n, 2, 0.01)
    @test_throws ArgumentError quiet() do
        computeCLV(flow_nl, flow_tan, flow_adj, zeros(n), () -> randn(n),
                   2.0, 0.25, 2, 2, 10)
    end
end
