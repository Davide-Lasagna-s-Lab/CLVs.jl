module CLVs

using Flows, 
      LinearAlgebra

export computeCLV, angles


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

    return times(mon_tan_cp) .- tspin, _process(samples(mon_adj_cp), samples(mon_tan_cp))
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

function _process(ψ⁺s::AbstractVector{Coupled{N, T}},
                  ψ⁻s::AbstractVector{Coupled{N, T}}) where {N, T}
    Γs = Coupled{N, T}[]

    for i in eachindex(ψ⁺s)
        # construct P matrix
        P = [dot(ψ⁺s[i][j], ψ⁻s[i][k]) for j = 1:N, k = 1:N]

        # compute LU decomposition, without pivoting
        A⁻ = inv(lu(P, Val(false)).U)

        # construct linear combination of backwards lyapunov vectors
        Γ = similar(ψ⁻s[1])
        for j in 1:N
            Γ[j] .*= 0
            for k in 1:N
                Γ[j] .+= ψ⁻s[i][k] .* A⁻[k, j]
            end
            Γ[j] ./= norm(Γ[j])
        end

        # store
        push!(Γs, Γ)
    end

    return Γs
end


""" Compute pairwise angles matrix between CLVs """
angles(Γ::Coupled{N}) where {N} =
    [acos.(abs.(clamp.(dot(Γ[i], Γ[j]), -1, 1))) for i = 1:N, j = 1:N]

dots(Γ::Coupled{N}) where {N} =
    [clamp.(dot(Γ[i], Γ[j]), -1, 1) for i = 1:N, j = 1:N]


end