module CLVs

using Flows, 
      LinearAlgebra

export computeCLV, angles


"""
    computeCLV(flow_nl, flow_tan, flow_adj, x, randx,
               tspin, T, Nsplits, N, oneevery=10) -> ts, Γs, Θs

Compute the first `N` covariant Lyapunov vectors (CLVs) and the adjoint CLVs
dual to them along the trajectory through `x`, over a window of length `T`
that follows a spin-up of length `tspin`. Return the sample times `ts`,
measured from the start of the window, the CLVs `Γs` and the adjoint CLVs
`Θs`. `Γs[i][j]` and `Θs[i][j]` are the `j`-th vectors at time `ts[i]`.

The arguments are:

  - `flow_nl`: flow of the nonlinear equations, called as `flow_nl(x, span)`
    and as `flow_nl(x, span, storage)`;
  - `flow_tan`, `flow_adj`: flows of the tangent and of the adjoint equations
    for `N` coupled vectors, along a stored trajectory, called as
    `flow_tan(ψ, storage, span, monitor)`. The adjoint one marches backward
    in time;
  - `x`: state at time zero, where the spin-up starts. It is modified;
  - `randx`: function returning a new random state at every call, for the
    initial tangent and adjoint vectors;
  - `tspin`: length of the spin-up of the tangent vectors, before the window,
    and of the adjoint vectors, after it;
  - `T`: length of the window, a whole number of sampling intervals;
  - `Nsplits`: number of chunks of the spin-up, whose trajectory is held in
    memory one chunk at a time;
  - `N`: number of vectors;
  - `oneevery`: number of time steps between two Gram-Schmidt
    orthonormalisations, and between two samples of the output.

The vectors are sorted by decreasing Lyapunov exponent. The CLVs have unit
`norm`, and the adjoint CLVs are scaled so that `dot(Θs[i][j], Γs[i][k])` is
one for `j == k` and zero otherwise. The adjoint equations must be the adjoint
of the tangent ones in the inner product `dot`. The time taken by each stage
of the calculation is printed to `stdout`.
"""
function computeCLV(flow_nl, flow_tan, flow_adj, x::X, randx, 
                    tspin::Real, T::Real, Nsplits::Int, N::Int, oneevery::Int=10) where {X}
    # construct the storage
    storage = RAMStorage(X);

    # initialise backward and forward lyapunov vectors,
    # to be marched forward and backward
    ψ⁻ = couple(ntuple(i->randx(), N)...)
    ψ⁺ = couple(ntuple(i->randx(), N)...)

    # this monitor will store the temporal variation of the Lyapunov 
    # exponents, plus perform the reorthogonalisation
    mon_tan = Monitor(ψ⁻, gs_ortho_Q!; oneevery=oneevery)
    mon_adj = Monitor(ψ⁺, gs_ortho_Q!; oneevery=oneevery)
    mon_tan_cp = Monitor(ψ⁻, gs_ortho_Q_cp!; oneevery=oneevery);
    mon_adj_cp = Monitor(ψ⁺, gs_ortho_Q_cp!; oneevery=oneevery);

    # define times
    tB = tspin
    tC = tB + T

    # initial time
    t = 0

    # length of the split
    Tchunk = tspin/Nsplits

    # store checkpoints here
    checkpoints = Tuple{Float64, X}[]

    ################ 
    # FORWARD 
    ################ 
    # spin up
    for i = 1:Nsplits
        @time flow_nl(x, (t, t+Tchunk), reset!(storage)); flush(stdout)
        @time flow_tan(ψ⁻, storage, (t, t+Tchunk), mon_tan); flush(stdout)
        t += Tchunk
    end

    # store initial condition of the point where the main calculations is done
    xB = copy(x)

    # main calculation done in a single go
    @time flow_nl(x, (tB, tC), reset!(storage)); flush(stdout)
    @time flow_tan(ψ⁻, storage, (tB, tC), mon_tan_cp); flush(stdout)

    # The tangent samples are counted from the start of the window and the
    # adjoint ones from its end. They fall on the same times only if the
    # window is a whole number of sampling intervals, that is if the tangent
    # samples are evenly spaced up to the last one.
    ts = times(mon_tan_cp)
    isapprox(ts[end] - ts[end-1], ts[2] - ts[1]; rtol=1e-6) ||
        throw(ArgumentError("T must be a whole number of sampling intervals, " *
                            "that is of `oneevery` time steps"))

    ################ 
    # BACKWARD 
    ################ 

    # update time
    t = tC

    # march NL equations forward and store checkpoints
    for i = 1:Nsplits
        # store checkpoint
        push!(checkpoints, (t, copy(x)))

        # march
        if i != Nsplits
            @time flow_nl(x, (t, t+Tchunk)); flush(stdout)
        end

        # advance time
        t += Tchunk
    end

    # backward integration - spin up
    for i = Nsplits:-1:1
        ti, xi = checkpoints[i]
        @time flow_nl(xi, (ti, ti+Tchunk), reset!(storage)); flush(stdout)
        @time flow_adj(ψ⁺, storage, (ti+Tchunk, ti), mon_adj); flush(stdout)
    end

    # main calculation done in a single go
    @time flow_nl(xB, (tB, tC), reset!(storage)); flush(stdout)
    @time flow_adj(ψ⁺, storage, (tC, tB), mon_adj_cp); flush(stdout)

    # reverse times for the adjoint case
    reverse!(times(mon_adj_cp))
    reverse!(samples(mon_adj_cp))

    # turn the samples of the two sets into the CLVs and the adjoint CLVs
    Γs, Θs = _process!(samples(mon_adj_cp), samples(mon_tan_cp))

    return times(mon_tan_cp) .- tspin, Γs, Θs
end

# re-orthogonalisation step for storage marching
function gs_ortho_Q!(t, ψs)
    N = length(ψs)
    R = zeros(N)
    for i = 1:N
        for j = 1:i-1
            r = dot(ψs[i], ψs[j])
            ψs[i] .-= r .* ψs[j]
        end
        R[i] = norm(ψs[i])
        ψs[i] ./= R[i]
    end
    return ψs
end

# copying version
gs_ortho_Q_cp!(t, ψs) = copy(gs_ortho_Q!(t, ψs))

# Turn the Gram-Schmidt vectors of the tangent equations, `ψ⁻s`, into the CLVs,
# and those of the adjoint equations, `ψ⁺s`, into the adjoint CLVs, in place,
# sample by sample. With P the matrix of the inner products between the two
# sets and P = LU its factorisation without pivoting, the CLVs are the columns
# of ψ⁻U⁻¹ and the adjoint CLVs those of ψ⁺L⁻ᵀ.
function _process!(ψ⁺s::AbstractVector{Coupled{N, T}},
                   ψ⁻s::AbstractVector{Coupled{N, T}}) where {N, T}
    for i in eachindex(ψ⁺s)
        ψ⁺ = ψ⁺s[i]
        ψ⁻ = ψ⁻s[i]

        # ---- inner products between the two sets ----
        P = [dot(ψ⁺[j], ψ⁻[k]) for j = 1:N, k = 1:N]

        # ---- coefficients of the two combinations ----
        # Both matrices are upper triangular, and B has unit diagonal.
        F = lu(P, NoPivot())
        A = inv(F.U)
        B = inv(F.L)'

        # ---- combine in place ----
        # Vector j needs the vectors from 1 to j only. Going from the last to
        # the first, each one is overwritten after its last use.
        for j in N:-1:1
            ψ⁻[j] .*= A[j, j]
            for k in 1:j-1
                ψ⁻[j] .+= ψ⁻[k] .* A[k, j]
                ψ⁺[j] .+= ψ⁺[k] .* B[k, j]
            end

            # CLV of unit norm, and adjoint CLV of unit inner product with it
            a = norm(ψ⁻[j])
            ψ⁻[j] ./= a
            ψ⁺[j] .*= a
        end
    end

    return ψ⁻s, ψ⁺s
end


""" Compute pairwise angles matrix between CLVs """
angles(Γ::Coupled{N}) where {N} =
    [acos.(abs.(clamp.(dot(Γ[i], Γ[j]), -1, 1))) for i = 1:N, j = 1:N]

dots(Γ::Coupled{N}) where {N} =
    [clamp.(dot(Γ[i], Γ[j]), -1, 1) for i = 1:N, j = 1:N]


end