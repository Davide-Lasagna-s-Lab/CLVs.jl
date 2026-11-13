# Makes the figures of the README in assets/ that show results of `computeCLV`, and prints every
# number that the README quotes: the example on the Lorenz equations, with the vectors on the
# attractor, their transversality and the components of a perturbation, and the convergence of
# the vectors with the spin-up on a linear system. The diagrams are in make_diagrams.jl, and the
# colours and theme of the figures in style.jl.
#
#     julia --project=examples -e 'using Pkg; Pkg.instantiate()'
#     julia --project=examples examples/make_figures.jl

using LinearAlgebra
using Printf
using Random
using Statistics
using Flows
using CLVs
using CairoMakie

include(joinpath(@__DIR__, "style.jl"))

const ASSETS = joinpath(@__DIR__, "..", "assets")

# `computeCLV` prints the time of each stage of the calculation. Run `f` without this output.
quiet(f) = redirect_stdout(f, devnull)

# The flows that `computeCLV` takes, for states of length n and N vectors: the flow of the
# nonlinear equations f, and those of the tangent and adjoint equations ftan and fadj along a
# stored trajectory, with steps of size Δt of the scheme `scheme`, RK4 or RK8.
function flows(f, ftan, fadj, n, N, Δt, scheme=RK4)
    z   = couple(ntuple(i -> zeros(n), N)...)
    dep = CallDependency(ntuple(i -> (i,), N)...)

    flow_nl  = flow(f, scheme(zeros(n)), TimeStepConstant(Δt))
    flow_tan = flow(couple(ntuple(i -> ftan, N)...), dep,
                    scheme(z, ContinuousMode(false)), TimeStepFromStorage(Δt))
    flow_adj = flow(couple(ntuple(i -> fadj, N)...), dep,
                    scheme(z, ContinuousMode(true)), TimeStepFromStorage(Δt))

    return flow_nl, flow_tan, flow_adj
end

# The trajectory over the window of `computeCLV`, from the state x at time zero, in a storage
# that holds it at every time step. The spin-up goes in the same chunks as in `computeCLV`, so
# that the trajectory is the same to the last bit.
function window(flow_nl, x, tspin, T, Nsplits; degree=3)
    u = copy(x)
    t = 0
    for i in 1:Nsplits
        flow_nl(u, (t, t + tspin/Nsplits))
        t += tspin/Nsplits
    end
    store = RAMStorage(u; degree=degree)
    flow_nl(u, (tspin, tspin + T), store)
    return store
end

# largest distance, over the samples, between the unit vectors vs and the directions of ws
function misalignment(vs, ws)
    worst = 0.0
    for (v, w) in zip(vs, ws)
        worst = max(worst, norm(v .- sign(dot(v, w)) .* w ./ norm(w)))
    end
    return worst
end


# ============================================================================ #
# The example of the README, on the Lorenz equations                           #
# ============================================================================ #

# Lorenz equations, with σ = 10, ρ = 28 and β = 8/3
function lorenz(t, x, dxdt)
    dxdt[1] = 10*(x[2] - x[1])
    dxdt[2] = 28*x[1] - x[2] - x[1]*x[3]
    dxdt[3] = -8/3*x[3] + x[1]*x[2]
    return dxdt
end

# tangent equations, the Jacobian at u applied to the tangent vector y
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

# number of vectors, time step, spin-up, window, chunks and sampling, as in the README
N        = 3
Δt       = 0.005
tspin    = 25.0
T        = 50.0
Nsplits  = 5
oneevery = 10

flow_nl, flow_tan, flow_adj = flows(lorenz, lorenz_tan, lorenz_adj, 3, N, Δt)

# a point on the attractor, kept because `computeCLV` modifies its argument
x  = flow_nl([1.0, 1.0, 20.0], (0, 20))
x₀ = copy(x)

Random.seed!(0)
ts, Γs, Θs, λs = quiet() do
    computeCLV(flow_nl, flow_tan, flow_adj, x, () -> randn(3), tspin, T, Nsplits, N, oneevery)
end

# ---- the numbers of the README ----
println("Example of the README\n")
println("julia> length(ts), ts[end]")
show(stdout, MIME"text/plain"(), (length(ts), ts[end]))
println("\n\njulia> [dot(Θs[1][j], Γs[1][k]) for j in 1:3, k in 1:3]")
show(stdout, MIME"text/plain"(), [dot(Θs[1][j], Γs[1][k]) for j in 1:3, k in 1:3])
println("\n\njulia> angles(Γs[1])")
show(stdout, MIME"text/plain"(), angles(Γs[1]))

c = [1/norm(Θ[j]) for Θ in Θs, j in 1:3]
println("\n\njulia> minimum(c; dims=1)")
show(stdout, MIME"text/plain"(), minimum(c; dims=1))
println("\n\njulia> maximum(c; dims=1)")
show(stdout, MIME"text/plain"(), maximum(c; dims=1))
println("\n\njulia> size(λs)")
show(stdout, MIME"text/plain"(), size(λs))
println("\n\njulia> mean(λs; dims=1)")
show(stdout, MIME"text/plain"(), mean(λs; dims=1))
println("\n\njulia> extrema(λs[:, 1])")
show(stdout, MIME"text/plain"(), extrema(λs[:, 1]))

v = [1.0, 2.0, 3.0]
a = [dot(Θs[1][j], v) for j in 1:3]
println("\n\njulia> norm(v .- sum(a[j] .* Γ[j] for j in 1:3))")
show(stdout, MIME"text/plain"(), norm(v .- sum(a[j] .* Γs[1][j] for j in 1:3)))
println("\n")

dual = maximum(abs(dot(Θ[j], Γ[k]) - (j == k)) for (Γ, Θ) in zip(Γs, Θs), j in 1:3, k in 1:3)
α    = [angles(Γ)[j, k] for Γ in Γs, (j, k) in ((1, 2), (1, 3), (2, 3))]
@printf "sampling interval %.3f, largest entry of ΘᵀΓ - I over the window %.1e\n" ts[2] - ts[1] dual
@printf "smallest angles between γ₁ γ₂, γ₁ γ₃ and γ₂ γ₃: %.4f %.4f %.4f\n" minimum(α; dims=1)...
@printf "sum of the mean growth rates %.4f, against the trace of the Jacobian %.4f\n" sum(mean(λs; dims=1)) -(10 + 1 + 8/3)
@printf "growth rates of the neutral CLV: from %.2f to %.2f\n" extrema(λs[:, 2])...

# ---- the trajectory over the window, at every time step and at the samples of the vectors ----
store = window(flow_nl, x₀, tspin, T, Nsplits)
ux = [s[1] for s in samples(store)]
uz = [s[3] for s in samples(store)]
uk = samples(store)[1:oneevery:end]
tk = times(store)[1:oneevery:end]
maximum(abs, tk .- tspin .- ts) < 1e-9 || error("the samples of the trajectory are misaligned")



# ============================================================================ #
# Figure: the vectors on the attractor                                         #
# ============================================================================ #

# The vectors are drawn over the first eight time units of the window, at one sample every three.
shown = findall(<=(8.0), ts)
picks = shown[1:3:end]
steps = 1:(last(shown) - 1)*oneevery + 1

# segments through the states of the trajectory along the projection on the plane (x, z) of the
# vectors j of the family V, of half length ℓ before the projection
function segments(V, j, ℓ)
    points = Point2f[]
    for k in picks
        w = ℓ .* V[k][j] ./ norm(V[k][j])
        push!(points, Point2f(uk[k][1] - w[1], uk[k][3] - w[3]), Point2f(uk[k][1] + w[1], uk[k][3] + w[3]))
    end
    return points
end

fig = Figure(size=(820, 500))
for (col, V, tag, names) in ((1, Γs, "(a)", [L"\gamma_1\ \mathrm{(unstable)}", L"\gamma_3\ \mathrm{(stable)}"]),
                             (2, Θs, "(b)", [L"\theta_1", L"\theta_3"]))
    ax = Axis(fig[1, col]; xlabel=L"x", ylabel=L"z", aspect=DataAspect(), limits=(-22, 22, 0, 52))
    lines!(ax, ux, uz; color=GRID, linewidth=1)
    lines!(ax, ux[steps], uz[steps]; color=AXIS, linewidth=2)

    # each segment on a ring of white
    for (j, c) in ((3, STABLE), (1, UNSTABLE))
        linesegments!(ax, segments(V, j, 2.0); color=:white, linewidth=5.5)
        linesegments!(ax, segments(V, j, 2.0); color=c, linewidth=2.4)
    end
    scatter!(ax, [uk[k][1] for k in picks], [uk[k][3] for k in picks]; color=INK, markersize=4.5)

    keys = [LineElement(color=UNSTABLE, linewidth=2.5), LineElement(color=STABLE, linewidth=2.5)]
    Legend(fig[0, col], keys, names; orientation=:horizontal, framevisible=false, labelcolor=INK2,
           tellwidth=false)
    panel!(ax, tag)
end
save(joinpath(ASSETS, "lorenz_vectors.png"), fig; px_per_unit=3)


# ============================================================================ #
# Figure: transversality along the trajectory                                  #
# ============================================================================ #

fig = Figure(size=(820, 400))
ax1 = Axis(fig[1, 1]; xlabel=L"t", ylabel=L"c_j", yscale=log10, yticks=LogTicks(-2:0))
# c₁ and c₂ nearly coincide. The second is drawn wider, under the first.
lines!(ax1, ts, c[:, 3]; color=STABLE, linewidth=1.5)
lines!(ax1, ts, c[:, 2]; color=NEUTRAL, linewidth=4.5, joinstyle=:round)
lines!(ax1, ts, c[:, 1]; color=UNSTABLE, linewidth=1.5)
xlims!(ax1, 0, 20)
ylims!(ax1, 0.015, 3)
Legend(fig[0, 1], [LineElement(color=UNSTABLE, linewidth=1.5), LineElement(color=NEUTRAL, linewidth=4.5),
                   LineElement(color=STABLE, linewidth=1.5)], [L"c_1", L"c_2", L"c_3"];
       orientation=:horizontal, framevisible=false, labelcolor=INK2, tellwidth=false)
panel!(ax1, "(a)")

# the length of the first adjoint CLV on the attractor, with the longest drawn last
ax2   = Axis(fig[1, 2]; xlabel=L"x", ylabel=L"z", aspect=DataAspect(), limits=(-22, 22, 0, 54))
order = sortperm(c[:, 1]; rev=true)
lines!(ax2, ux, uz; color=GRID, linewidth=1)
sc = scatter!(ax2, [s[1] for s in uk][order], [s[3] for s in uk][order]; color=1 ./ c[order, 1],
              colormap=SEQUENTIAL, colorscale=log10, colorrange=(1, 50), markersize=6)
Colorbar(fig[1, 3], sc; label=L"\Vert \theta_1 \Vert = 1/c_1", ticks=[1, 2, 5, 10, 20, 50],
         labelcolor=INK, ticklabelcolor=INK2, tickcolor=AXIS, spinewidth=0)
panel!(ax2, "(b)")
colsize!(fig.layout, 1, Relative(0.56))
save(joinpath(ASSETS, "lorenz_transversality.png"), fig; px_per_unit=3)

near = argmin(c[:, 1])
@printf "c₁ is smallest at t = %.2f, at the state (%.2f, %.2f, %.2f)\n" ts[near] uk[near]...
@printf "fraction of the samples with c₁ below 0.1: %.3f\n" count(<(0.1), c[:, 1])/length(ts)


# ============================================================================ #
# Figure: the components of a perturbation on the CLVs                         #
# ============================================================================ #

# A perturbation is marched with the tangent flow along the window, and its components on the
# CLVs are its inner products with the adjoint CLVs. The three copies are a set of N vectors.
Random.seed!(1)
y₀ = randn(3)
y  = couple(copy(y₀), copy(y₀), copy(y₀))
b  = zeros(length(ts), 3)
for i in eachindex(ts)
    i > 1 && flow_tan(y, store, (tk[i-1], tk[i]))
    for j in 1:3
        b[i, j] = abs(dot(Θs[i][j], y[1]))
    end
end

# growth of each component from its value at the start of the window
g  = b ./ b[1:1, :]
λs = (0.906, 0.0, -14.57)

fig = Figure(size=(820, 360))
ax1 = Axis(fig[1, 1]; xlabel=L"t",
           ylabel=L"|\langle \theta_j(t), y(t) \rangle|\ /\ |\langle \theta_j(0), y(0) \rangle|",
           yscale=log10, yticks=LogTicks(-2:2:8))
lines!(ax1, ts, exp.(λs[1] .* ts); color=MUTED, linewidth=1)
lines!(ax1, ts, g[:, 1]; color=UNSTABLE, label=L"j = 1")
lines!(ax1, ts, g[:, 2]; color=NEUTRAL, label=L"j = 2")
text!(ax1, 11.2, exp(λs[1]*11.2)*8; text=L"e^{0.906 t}", align=(:right, :bottom), color=INK2, fontsize=15)
xlims!(ax1, 0, 16)
ylims!(ax1, 1e-2, 1e8)
axislegend(ax1; position=:lt, margin=(34, 10, 10, 10))
panel!(ax1, "(a)")

ax2 = Axis(fig[1, 2]; xlabel=L"t",
           ylabel=L"|\langle \theta_3(t), y(t) \rangle|\ /\ |\langle \theta_3(0), y(0) \rangle|",
           yscale=log10, yticks=LogTicks(-12:3:0))
lines!(ax2, ts, exp.(λs[3] .* ts); color=MUTED, linewidth=1)
lines!(ax2, ts, g[:, 3]; color=STABLE, label=L"j = 3")
text!(ax2, 1.02, exp(λs[3]*1.0)*4; text=L"e^{-14.57 t}", align=(:left, :bottom), color=INK2, fontsize=15)
xlims!(ax2, 0, 1.8)
ylims!(ax2, 1e-12, 1e3)
axislegend(ax2; position=:rt)
panel!(ax2, "(b)")
colsize!(fig.layout, 1, Relative(0.58))
save(joinpath(ASSETS, "lorenz_projection.png"), fig; px_per_unit=3)

i₁ = findfirst(>=(16), ts)
i₃ = findfirst(>=(2), ts)
@printf "growth rates of the components: j = 1 over 16 time units %.3f, j = 3 over 2 time units %.2f\n" log(g[i₁, 1])/ts[i₁] log(g[i₃, 3])/ts[i₃]
@printf "the component j = 2 stays between %.2f and %.2f over the first 16 time units\n" extrema(g[1:i₁, 2])...


# ============================================================================ #
# Accuracy: the time step and the number of vectors                            #
# ============================================================================ #

println("\nAccuracy")

# the neutral CLV against the direction of the flow, with the same sampling, for RK4 with three
# time steps and for RK8 with two degrees of the interpolation of the stored trajectory
for (name, scheme, h, degree) in (("RK4", RK4, 0.005, 3), ("RK4", RK4, 0.0025, 3), ("RK4", RK4, 0.00125, 3),
                                  ("RK8", RK8, 0.005, 3), ("RK8", RK8, 0.005, 7))
    every      = round(Int, 0.05/h)
    F, Ft, Fa  = flows(lorenz, lorenz_tan, lorenz_adj, 3, 3, h, scheme)
    w          = F([1.0, 1.0, 20.0], (0, 20))
    states     = samples(window(F, w, tspin, T, Nsplits; degree=degree))[1:every:end]
    Random.seed!(0)
    _, G, _ = quiet() do
        computeCLV(F, Ft, Fa, w, () -> randn(3), tspin, T, Nsplits, 3, every; degree=degree)
    end
    off = misalignment([Γ[2] for Γ in G], [lorenz(0.0, s, zeros(3)) for s in states])
    @printf "%s, time step %.5f, degree %d: neutral CLV against the direction of the flow %.1e\n" name h degree off
end

# the first two vectors of each family, with two vectors instead of three
F, Ft, Fa = flows(lorenz, lorenz_tan, lorenz_adj, 3, 2, Δt)
Random.seed!(3)
_, Γ₂, Θ₂ = quiet() do
    computeCLV(F, Ft, Fa, copy(x₀), () -> randn(3), tspin, T, Nsplits, 2, oneevery)
end
errΓ = 0.0
errΘ = 0.0
for i in eachindex(ts), j in 1:2
    s = sign(dot(Γ₂[i][j], Γs[i][j]))
    global errΓ = max(errΓ, norm(Γ₂[i][j] .- s .* Γs[i][j]))
    global errΘ = max(errΘ, norm(Θ₂[i][j] .- s .* Θs[i][j]) / norm(Θs[i][j]))
end
@printf "two vectors against three: CLVs %.1e, adjoint CLVs %.1e\n" errΓ errΘ


# ============================================================================ #
# Figure: convergence with the spin-up, on a linear system                     #
# ============================================================================ #

# The fixed point x = 0 of x' = A x, where the CLVs are the right eigenvectors of A and the
# adjoint CLVs the left ones. The gaps between the four exponents are 1.2, 0.7 and 0.35.
n = 4
λ = [1.0, -0.2, -0.9, -1.25]
V = Matrix(1.0I, n, n) .+ 0.4 .* [sin(i + 2j) for i in 1:n, j in 1:n]
V = V ./ sqrt.(sum(abs2, V; dims=1))
W = Matrix(inv(V)')
A = V * Diagonal(λ) / V

linear(t, x, dxdt)        = mul!(dxdt, A, x)
linear_tan(t, u, y, dydt) = mul!(dydt, A, y)
linear_adj(t, u, w, dwdt) = mul!(dwdt, A', w)

flow_nl, flow_tan, flow_adj = flows(linear, linear_tan, linear_adj, n, 3, 0.01)

# largest error of the first three vectors of each family over a window of one time unit, with
# the same random initial vectors for every length of the spin-up
spins = 4.0:4.0:60.0
errΓ  = zeros(length(spins), 3)
errΘ  = zeros(length(spins), 3)
for (k, τ) in enumerate(spins)
    Random.seed!(2)
    _, G, H = quiet() do
        computeCLV(flow_nl, flow_tan, flow_adj, zeros(n), () -> randn(n), τ, 1.0, 4, 3, 10)
    end
    for (Γ, Θ) in zip(G, H), j in 1:3
        s          = sign(dot(Γ[j], V[:, j]))
        errΓ[k, j] = max(errΓ[k, j], norm(Γ[j] .- s .* V[:, j]))
        errΘ[k, j] = max(errΘ[k, j], norm(Θ[j] .- s .* W[:, j]) / norm(W[:, j]))
    end
end

# the vector j converges at the smaller of the gaps on its two sides
rates = (1.2, 0.7, 0.35)

fig = Figure(size=(820, 380))
for (col, err, tag, name) in ((1, errΓ, "(a)", L"\mathrm{error\ of\ } \gamma_j"),
                              (2, errΘ, "(b)", L"\mathrm{error\ of\ } \theta_j"))
    ax = Axis(fig[1, col]; xlabel=L"t_\mathrm{spin}", ylabel=name, yscale=log10, yticks=LogTicks(-16:4:0))
    for j in 1:3
        # the predicted decay, through the measured error at the second spin-up
        lines!(ax, spins, err[2, j] .* exp.(-rates[j] .* (spins .- spins[2])); color=MUTED, linewidth=1)
        scatterlines!(ax, spins, max.(err[:, j], 1e-17); color=SLOTS[j], label=L"j = %$j", markersize=8,
                      strokecolor=:white, strokewidth=1.5)

        # the rate, beside the line: below it for the first vector, above it for the others
        τ     = (21, 33, 44)[j]
        shift = (1/40, 40, 40)[j]
        align = j == 1 ? (:right, :top) : (:left, :bottom)
        text!(ax, τ, shift*err[2, j]*exp(-rates[j]*(τ - spins[2]));
              text=L"e^{-%$(rates[j]) t_\mathrm{spin}}", align, color=INK2, fontsize=14)
    end
    ylims!(ax, 1e-16, 1e3)
    xlims!(ax, 0, 63)
    col == 1 && Legend(fig[0, 1:2], ax; orientation=:horizontal, framevisible=false, labelcolor=INK2,
                       tellwidth=false)
    panel!(ax, tag)
end
save(joinpath(ASSETS, "spinup.png"), fig; px_per_unit=3)

println("\nConvergence with the spin-up, errors of γ₁ γ₂ γ₃ and of θ₁ θ₂ θ₃")
for (k, τ) in enumerate(spins)
    @printf "tspin %4.0f   %.1e %.1e %.1e   %.1e %.1e %.1e\n" τ errΓ[k, :]... errΘ[k, :]...
end
for (name, err) in (("γ", errΓ), ("θ", errΘ)), j in 1:3
    # least squares fit of the decay, over the lengths at which the error is above rounding
    k    = findall(>(1e-12), err[:, j])
    fit  = -([ones(length(k)) spins[k]] \ log.(err[k, j]))[2]
    @printf "measured rate of %s%d up to tspin = %.0f: %.2f, predicted %.2f\n" name j spins[k[end]] fit rates[j]
end

println("\nwritten assets/lorenz_vectors.png, assets/lorenz_transversality.png, assets/lorenz_projection.png and assets/spinup.png")
