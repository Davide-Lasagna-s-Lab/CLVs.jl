# Draws the diagrams of the README in assets/: the covariance of the vectors, the dual basis, the
# transversality, the convergence of the directions, the nested subspaces, the factorisation, the
# algorithm and the sampling. They are drawings with exact geometry and involve no integration.
# The figures with the results of `computeCLV` are in make_figures.jl.
#
#     julia --project=examples examples/make_diagrams.jl
#
# With the names of some diagrams as arguments, only these are drawn, for instance
#
#     julia --project=examples examples/make_diagrams.jl duality subspaces
#
# In every diagram the first, second and third vector of a family are orange, green and blue.
# The vectors of the tangent equations, the CLVs γ and the Gram-Schmidt vectors ψ⁻, are drawn
# with a solid line, and those of the adjoint equations, θ and ψ⁺, with a dashed line.

using LinearAlgebra
using CairoMakie

include(joinpath(@__DIR__, "style.jl"))

const ASSETS = joinpath(@__DIR__, "..", "assets")

# the dash pattern of the adjoint vectors
const DASH = Linestyle([0.0, 2.2, 3.6])


# ============================================================================ #
# Drawing helpers                                                              #
# ============================================================================ #

# unit vector at an angle in degrees, and the vector a quarter turn ahead of d
unitvec(a) = Point2f(cosd(a), sind(a))
turn(d)    = Point2f(-d[2], d[1])

# A figure of w by h units with one axis that covers it, in which a length of one is one unit of
# the figure. The diagrams are drawn in these units.
function canvas(w, h)
    fig = Figure(size=(w, h), figure_padding=0)
    ax  = Axis(fig[1, 1]; limits=(0, w, 0, h), aspect=DataAspect())
    hidedecorations!(ax)
    hidespines!(ax)
    return fig, ax
end

# an arrow from p to q, with a filled head
function arrow!(ax, p, q; color=INK, linewidth=2.6, linestyle=:solid, head=12)
    len  = norm(q - p)
    d    = (q - p) / len
    head = min(head, 0.5*len)
    base = q - head*d

    lines!(ax, Point2f[p, base + 0.3*head*d]; color, linewidth, linestyle, linecap=:butt)
    poly!(ax, Point2f[q, base + 0.38*head*turn(d), base - 0.38*head*turn(d)]; color, strokewidth=0)
    return nothing
end

# a segment from p to q
segment!(ax, p, q; color=MUTED, linewidth=1, linestyle=:solid) =
    lines!(ax, Point2f[p, q]; color, linewidth, linestyle)

# the mark of a right angle at p, between the unit vectors d₁ and d₂
rightangle!(ax, p, d₁, d₂; size=9, color=MUTED) =
    lines!(ax, Point2f[p + size*d₁, p + size*(d₁ + d₂), p + size*d₂]; color, linewidth=1)

# an arc of radius r about c, from the angle a₀ to the angle a₁ in degrees
arc!(ax, c, r, a₀, a₁; color=MUTED, linewidth=1, linestyle=:solid) =
    lines!(ax, Point2f[c + r*unitvec(a) for a in range(a₀, a₁; length=90)]; color, linewidth, linestyle)

# a label in mathematics or in plain text, centred at p unless told otherwise
label!(ax, p, s; align=(:center, :center), fontsize=18, color=INK, font=:regular) =
    text!(ax, p; text=s, align, fontsize, color, font)

# the label of a panel, with its top left corner at p
panellabel!(ax, p, s) = label!(ax, p, s; align=(:left, :top), fontsize=15, font=:bold)


# ============================================================================ #
# Covariance: the tangent flow maps each CLV onto itself                       #
# ============================================================================ #

function covariance()
    fig, ax = canvas(820, 305)

    # ---- the trajectory, a cubic arch, and its direction ----
    b₀ = Point2f(20, 35)
    b₁ = Point2f(250, 237)
    b₂ = Point2f(560, 227)
    b₃ = Point2f(800, 35)

    curve(τ)   = (1 - τ)^3*b₀ + 3*(1 - τ)^2*τ*b₁ + 3*(1 - τ)*τ^2*b₂ + τ^3*b₃
    slope(τ)   = 3*(1 - τ)^2*(b₁ - b₀) + 6*(1 - τ)*τ*(b₂ - b₁) + 3*τ^2*(b₃ - b₂)
    heading(τ) = atand(slope(τ)[2], slope(τ)[1])

    lines!(ax, Point2f[curve(τ) for τ in range(0, 0.985; length=200)]; color=AXIS, linewidth=2.5)
    arrow!(ax, curve(0.97), curve(1.0); color=AXIS, linewidth=2.5, head=14)

    # ---- the states at the two times ----
    τₛ = 0.2
    τₜ = 0.8
    uₛ = curve(τₛ)
    uₜ = curve(τₜ)
    aₛ = heading(τₛ)
    aₜ = heading(τₜ)

    # ---- unit CLVs at time s, and their images at time t ----
    # The neutral vector is along the trajectory. The angles are measured from it.
    len = 66
    γₛ  = (unitvec(aₛ + 72), unitvec(aₛ), unitvec(aₛ + 141))
    γₜ  = (unitvec(aₜ + 58), unitvec(aₜ), unitvec(aₜ + 152))
    μ   = (1.75, 1.12, 0.46)

    # a cell of perturbations between the first and the third vector, and its image
    cellₛ = Point2f[uₛ, uₛ + 0.5*len*γₛ[1], uₛ + 0.5*len*(γₛ[1] + γₛ[3]), uₛ + 0.5*len*γₛ[3]]
    cellₜ = Point2f[uₜ, uₜ + 0.5*len*μ[1]*γₜ[1], uₜ + 0.5*len*(μ[1]*γₜ[1] + μ[3]*γₜ[3]),
             uₜ + 0.5*len*μ[3]*γₜ[3]]
    poly!(ax, cellₛ; color=(MUTED, 0.22), strokewidth=0)
    poly!(ax, cellₜ; color=(MUTED, 0.22), strokewidth=0)

    for j in (3, 2, 1)
        arrow!(ax, uₛ, uₛ + len*γₛ[j]; color=SLOTS[j])
        arrow!(ax, uₜ, uₜ + len*μ[j]*γₜ[j]; color=SLOTS[j])
    end
    scatter!(ax, Point2f[uₛ, uₜ]; color=INK, markersize=9, strokecolor=:white, strokewidth=2)

    # ---- labels ----
    label!(ax, uₛ + Point2f(8, -20), L"u(s)")
    label!(ax, uₜ + Point2f(-8, -20), L"u(t)")

    label!(ax, uₛ + (len + 18)*γₛ[1] + Point2f(8, 0), L"\gamma_1(s)")
    label!(ax, uₛ + len*γₛ[2] + Point2f(28, -17), L"\gamma_2(s)")
    label!(ax, uₛ + (len + 14)*γₛ[3] + Point2f(-26, 2), L"\gamma_3(s)")

    label!(ax, uₜ + len*μ[1]*γₜ[1] + Point2f(-6, 20), L"\mu_1 \gamma_1(t)")
    label!(ax, uₜ + len*μ[2]*γₜ[2] + Point2f(52, 14), L"\mu_2 \gamma_2(t)")
    label!(ax, uₜ + len*μ[3]*γₜ[3] + Point2f(-32, 38), L"\mu_3 \gamma_3(t)")

    # ---- the propagator, from the vectors at s to those at t ----
    c  = Point2f(410, -175)
    r  = 452
    a₀ = 113
    a₁ = 70.5
    arc!(ax, c, r, a₀, a₁ + 1.5; color=INK2, linewidth=1.3)
    arrow!(ax, c + r*unitvec(a₁ + 1.6), c + r*unitvec(a₁); color=INK2, linewidth=1.3, head=11)
    label!(ax, Point2f(410, 291), L"Y(t, s)")

    save(joinpath(ASSETS, "covariance.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Duality: the adjoint CLVs are the basis dual to the CLVs                     #
# ============================================================================ #

function duality()
    fig, ax = canvas(820, 350)

    # ---- two CLVs of the plane and their adjoints ----
    a₁ = 8.0                  # angle of γ₁ in degrees
    a₂ = 62.0                 # angle of γ₂
    γ₁ = unitvec(a₁)
    γ₂ = unitvec(a₂)
    θ₁ = unitvec(a₂ - 90) / sind(a₂ - a₁)
    θ₂ = unitvec(a₁ + 90) / sind(a₂ - a₁)

    # ---- panel (a): the two bases ----
    o = Point2f(150, 140)     # origin
    U = 112                   # length of a unit vector

    # the lines spanned by the CLVs
    segment!(ax, o - 0.55*U*γ₁, o + 1.75*U*γ₁; color=AXIS)
    segment!(ax, o - 0.55*U*γ₂, o + 1.6*U*γ₂; color=AXIS)

    # each adjoint is normal to the line of the other CLV
    rightangle!(ax, o, -γ₂, θ₁/norm(θ₁))
    rightangle!(ax, o, -γ₁, θ₂/norm(θ₂))

    # the foot of the normal from the tip of each adjoint to its CLV is the tip of the CLV
    segment!(ax, o + U*θ₁, o + U*γ₁; color=MUTED, linestyle=:dot, linewidth=1.5)
    segment!(ax, o + U*θ₂, o + U*γ₂; color=MUTED, linestyle=:dot, linewidth=1.5)
    rightangle!(ax, o + U*γ₁, γ₁, -turn(γ₁); size=8)
    rightangle!(ax, o + U*γ₂, γ₂, turn(γ₂); size=8)

    arrow!(ax, o, o + U*θ₁; color=UNSTABLE, linestyle=DASH)
    arrow!(ax, o, o + U*θ₂; color=STABLE, linestyle=DASH)
    arrow!(ax, o, o + U*γ₁; color=UNSTABLE)
    arrow!(ax, o, o + U*γ₂; color=STABLE)

    label!(ax, o + 0.62*U*γ₁ + Point2f(2, -15), L"\gamma_1")
    label!(ax, o + 0.62*U*γ₂ + Point2f(17, -6), L"\gamma_2")
    label!(ax, o + U*θ₁ + Point2f(16, -8), L"\theta_1")
    label!(ax, o + U*θ₂ + Point2f(-16, 8), L"\theta_2")
    label!(ax, Point2f(292, 300), L"\langle \theta_i, \gamma_j \rangle = \delta_{ij}")
    panellabel!(ax, Point2f(16, 336), "(a)")

    # ---- panel (b): the components of a perturbation ----
    o = Point2f(455, 62)
    U = 74
    c = (2.25, 1.55)          # components of y on the two CLVs
    y = c[1]*γ₁ + c[2]*γ₂

    # the lines of constant ⟨θ₁, y⟩ are parallel to γ₂, and those of constant ⟨θ₂, y⟩ to γ₁
    for k in 0:3
        segment!(ax, o + U*(k*γ₁ - 0.12*γ₂), o + U*(k*γ₁ + 2.75*γ₂); color=GRID)
        label!(ax, o + U*(k*γ₁ + 2.75*γ₂) + Point2f(5, 12), string(k); fontsize=13, color=INK2)
    end
    for k in 0:2
        segment!(ax, o + U*(k*γ₂ - 0.12*γ₁), o + U*(k*γ₂ + 3.45*γ₁); color=GRID)
        label!(ax, o + U*(k*γ₂ + 3.45*γ₁) + Point2f(14, 1), string(k); fontsize=13, color=INK2)
    end
    label!(ax, o + U*(1.5*γ₁ + 2.75*γ₂) + Point2f(5, 36), L"\langle \theta_1, y \rangle")
    label!(ax, o + U*(1.0*γ₂ + 3.45*γ₁) + Point2f(42, 28), L"\langle \theta_2, y \rangle")

    # the perturbation and its two components
    segment!(ax, o + U*c[1]*γ₁, o + U*y; color=MUTED, linestyle=:dot, linewidth=1.5)
    segment!(ax, o + U*c[2]*γ₂, o + U*y; color=MUTED, linestyle=:dot, linewidth=1.5)
    arrow!(ax, o, o + U*c[1]*γ₁; color=UNSTABLE)
    arrow!(ax, o, o + U*c[2]*γ₂; color=STABLE)
    arrow!(ax, o, o + U*y; color=INK)

    label!(ax, o + U*y + Point2f(12, 10), L"y")
    label!(ax, o + 0.5*U*c[1]*γ₁ + Point2f(8, -19), L"\langle \theta_1, y \rangle \gamma_1")
    label!(ax, o + 0.62*U*c[2]*γ₂ + Point2f(-56, 4), L"\langle \theta_2, y \rangle \gamma_2")
    panellabel!(ax, Point2f(396, 336), "(b)")

    save(joinpath(ASSETS, "duality.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Transversality: the length of the adjoint CLV and the angle to the others    #
# ============================================================================ #

function transversality()
    fig, ax = canvas(820, 330)

    U = 68                    # length of a unit vector
    y = 82                    # height of the span of the other CLVs

    for (x, φ) in ((120, 75), (345, 45), (600, 19))
        o  = Point2f(x, y)
        γ₁ = unitvec(φ)
        γ₂ = unitvec(0)
        θ₁ = unitvec(90) / sind(φ)

        # the span of the other CLVs
        segment!(ax, o - Point2f(62, 0), o + Point2f(150, 0); color=AXIS, linewidth=1.5)

        # the distance of the tip of γ₁ from this span is c₁
        foot = Point2f(o[1] + U*cosd(φ), y)
        segment!(ax, foot, o + U*γ₁; color=INK2, linewidth=1.2)
        rightangle!(ax, foot, unitvec(180), unitvec(90); size=7)
        label!(ax, foot + Point2f(13, 0.5*U*sind(φ)), L"c_1"; fontsize=16)

        # the foot of the normal from the tip of θ₁ to γ₁ is the tip of γ₁
        segment!(ax, o + U*θ₁, o + U*γ₁; color=MUTED, linestyle=:dot, linewidth=1.5)

        # the angle between γ₁ and the span
        arc!(ax, o, 26, 0, φ; color=INK2)

        arrow!(ax, o, o + U*θ₁; color=UNSTABLE, linestyle=DASH)
        arrow!(ax, o, o + U*γ₁; color=UNSTABLE)
        arrow!(ax, o, o + U*γ₂; color=STABLE)

        label!(ax, o + U*θ₁ + Point2f(-17, -4), L"\theta_1")
        label!(ax, o + U*γ₁ + Point2f(15, 9), L"\gamma_1")
        label!(ax, o + U*γ₂ + Point2f(8, -15), L"\gamma_2")

        c = round(sind(φ); digits=2)
        label!(ax, o + Point2f(44, -52), L"c_1 = %$c"; fontsize=16)
    end

    save(joinpath(ASSETS, "transversality.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Directions: tangent vectors turn to γ₁, adjoint vectors to θ₁                #
# ============================================================================ #

function directions()
    fig, ax = canvas(820, 320)

    # ---- a linear system of the plane, with exponents 1 and -1 ----
    a₁ = 20.0                 # angle of γ₁ in degrees
    a₂ = 65.0                 # angle of γ₂
    V  = [cosd(a₁) cosd(a₂); sind(a₁) sind(a₂)]
    W  = inv(V)'              # the adjoint CLVs, in its columns
    b₁ = atand(W[2, 1], W[1, 1])
    b₂ = atand(W[2, 2], W[1, 2])

    # direction after a time τ of the tangent vector, or of the adjoint vector marched backward,
    # that starts at the angle a
    tangent(a, τ) = normalize(V * Diagonal([exp(τ), exp(-τ)]) * (V \ [cosd(a), sind(a)]))
    adjoint(a, τ) = normalize(W * Diagonal([exp(τ), exp(-τ)]) * (W \ [cosd(a), sind(a)]))

    R  = 108                  # radius of the circle of the directions
    τs = (0.0, 0.45, 0.9, 1.7)

    for (o, attr, rep, march, a₀, style, name) in
            ((Point2f(215, 178), a₁, a₂, tangent, 56.0, :solid, "tangent vectors, forward in time"),
             (Point2f(605, 178), b₁, b₂, adjoint, -61.0, DASH, "adjoint vectors, backward in time"))

        # ---- the circle of the directions, and the way they turn ----
        arc!(ax, o, R, 0, 360; color=AXIS, linewidth=1.5)
        for (r, a) in ((rep, attr), (rep, attr + 180), (rep + 180, attr + 180), (rep + 180, attr + 360))
            for f in (0.38, 0.68)
                α = r + f*(a - r)
                p = o + R*unitvec(α)
                d = sign(a - r) * turn(unitvec(α))
                poly!(ax, Point2f[p + 7*d, p - 5*d + 5*turn(d), p - 5*d - 5*turn(d)];
                      color=MUTED, strokewidth=0)
            end
        end

        # ---- the two directions that do not turn ----
        segment!(ax, o - R*unitvec(attr), o; color=(UNSTABLE, 0.45), linewidth=1.5, linestyle=style)
        segment!(ax, o - R*unitvec(rep), o; color=(STABLE, 0.45), linewidth=1.5, linestyle=style)

        # ---- one vector at four times, from light to dark ----
        for (τ, shade) in zip(τs, (0.25, 0.4, 0.6, 0.9))
            arrow!(ax, o, o + R*Point2f(march(a₀, τ)); color=(INK2, shade), linewidth=2, head=10,
                   linestyle=style)
        end

        arrow!(ax, o, o + R*unitvec(rep); color=STABLE, linestyle=style)
        arrow!(ax, o, o + R*unitvec(attr); color=UNSTABLE, linestyle=style)
        scatter!(ax, Point2f[o + R*unitvec(attr), o - R*unitvec(attr)]; color=UNSTABLE, markersize=11,
                 strokecolor=:white, strokewidth=2)
        scatter!(ax, Point2f[o + R*unitvec(rep), o - R*unitvec(rep)]; color=STABLE, markersize=11,
                 strokecolor=:white, strokewidth=2)

        label!(ax, o + Point2f(0, -R - 34), name; fontsize=14, color=INK2)
    end

    # ---- labels of panel (a) ----
    o = Point2f(215, 178)
    label!(ax, o + (R + 20)*unitvec(a₁), L"\gamma_1")
    label!(ax, o + (R + 18)*unitvec(a₂), L"\gamma_2")
    label!(ax, o + (R + 22)*Point2f(tangent(56.0, 0.5)), L"\psi_1^-")
    panellabel!(ax, Point2f(16, 306), "(a)")

    # ---- labels of panel (b), with the CLV γ₂ to which θ₁ is normal ----
    o = Point2f(605, 178)
    arrow!(ax, o, o + R*unitvec(a₂); color=(STABLE, 0.5), linewidth=1.8, head=10)
    rightangle!(ax, o, unitvec(a₂), unitvec(b₁))
    label!(ax, o + (R + 18)*unitvec(a₂), L"\gamma_2")
    label!(ax, o + (R + 20)*unitvec(b₁), L"\theta_1")
    label!(ax, o + (R + 18)*unitvec(b₂), L"\theta_2")
    label!(ax, o + (R + 22)*Point2f(adjoint(-61.0, 0.5)), L"\psi_1^+")
    panellabel!(ax, Point2f(406, 306), "(b)")

    save(joinpath(ASSETS, "directions.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Subspaces: each vector is the intersection of two subspaces                  #
# ============================================================================ #

function subspaces()
    fig, ax = canvas(820, 350)

    # ---- three CLVs of the space and their adjoints ----
    γ₁ = [1.0, 0.0, 0.0]
    γ₂ = [cosd(50), sind(50), 0.0]
    m  = [-sind(50), cosd(50), 0.0]
    q  = normalize([0.0, 0.0, 1.0] + tand(15)*m)
    γ₃ = normalize(cosd(20)*q + sind(20)*γ₂)
    Γ  = hcat(γ₁, γ₂, γ₃)
    Θ  = inv(Γ)'              # the adjoint CLVs, in its columns

    # ---- the view, from the azimuth α and the elevation ε in degrees ----
    α = 0
    ε = 45
    right = [-sind(α), cosd(α), 0.0]
    up    = [-sind(ε)*cosd(α), -sind(ε)*sind(α), cosd(ε)]
    front = [cosd(ε)*cosd(α), cosd(ε)*sind(α), sind(ε)]

    # ---- one panel ----
    # The vectors v₁, v₂ and v₃ of one family, and the first Gram-Schmidt vector n of the other
    # set. The plane of v₁ and v₂ and the plane normal to n, which holds v₂ and v₃, meet along v₂.
    function draw(o, U, v₁, v₂, v₃, n, style, other, names, offsets)
        project(v) = o + U*Point2f(dot(v, right), dot(v, up))

        # directions in the two planes, normal to their intersection
        e₂ = normalize(v₂)
        e₁ = normalize(v₁ - dot(v₁, e₂)*e₂)
        e₃ = normalize(v₃ - dot(v₃, e₂)*e₂)

        # the two halves of each plane, from the farthest to the nearest
        halves = [(e₁, -0.6, 0.0, style), (e₁, 0.0, 1.25, style),
                  (e₃, -0.5, 0.0, other), (e₃, 0.0, 1.3, other)]
        sort!(halves; by = h -> dot((h[2] + h[3])/2 * h[1], front))
        for (e, t₀, t₁, ls) in halves
            corners = Point2f[project(-0.55*e₂ + t₀*e), project(1.5*e₂ + t₀*e),
                              project(1.5*e₂ + t₁*e), project(-0.55*e₂ + t₁*e)]
            poly!(ax, corners; color=(MUTED, 0.16), strokewidth=0)
            lines!(ax, corners[[2, 3, 4, 1]]; color=MUTED, linewidth=1, linestyle=ls)
        end
        segment!(ax, project(-0.55*e₂), project(1.5*e₂); color=(NEUTRAL, 0.5), linewidth=1.5)

        # n is normal to the second plane
        s = 0.13
        lines!(ax, Point2f[project(s*normalize(n)), project(s*(normalize(n) + e₃)), project(s*e₃)];
               color=MUTED, linewidth=1)

        arrow!(ax, project(zeros(3)), project(normalize(n)); color=INK, linestyle=other)
        arrow!(ax, project(zeros(3)), project(v₃); color=STABLE, linestyle=style)
        arrow!(ax, project(zeros(3)), project(v₁); color=UNSTABLE, linestyle=style)
        arrow!(ax, project(zeros(3)), project(v₂); color=NEUTRAL, linestyle=style)

        for (v, name, offset) in zip((v₁, v₂, v₃, normalize(n)), names, offsets)
            label!(ax, project(v) + offset, name)
        end
        return nothing
    end

    # ---- panel (a): the CLV γ₂ ----
    o = Point2f(200, 172)
    draw(o, 112, γ₁, γ₂, γ₃, Θ[:, 1], :solid, DASH,
         (L"\gamma_1 = \psi_1^-", L"\gamma_2", L"\gamma_3", L"\psi_1^+"),
         (Point2f(40, -10), Point2f(18, 4), Point2f(-16, 10), Point2f(-8, -20)))
    label!(ax, o + Point2f(-95, -112), L"\mathrm{span}(\psi_1^-, \psi_2^-)"; fontsize=16)
    label!(ax, o + Point2f(128, 117), L"\mathrm{normal\ to\ } \psi_1^+"; fontsize=16)
    panellabel!(ax, Point2f(16, 336), "(a)")

    # ---- panel (b): the adjoint CLV θ₂ ----
    o = Point2f(590, 172)
    draw(o, 112, Θ[:, 1], Θ[:, 2], Θ[:, 3], γ₁, DASH, :solid,
         (L"\theta_1", L"\theta_2", L"\theta_3", L"\psi_1^- = \gamma_1"),
         (Point2f(-16, -8), Point2f(18, 4), Point2f(-16, 10), Point2f(54, -4)))
    label!(ax, o + Point2f(-72, -124), L"\mathrm{span}(\psi_1^+, \psi_2^+)"; fontsize=16)
    label!(ax, o + Point2f(146, 92), L"\mathrm{normal\ to\ } \psi_1^-"; fontsize=16)
    panellabel!(ax, Point2f(416, 336), "(b)")

    save(joinpath(ASSETS, "subspaces.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Factorisation: the shapes of the matrices of the recombination               #
# ============================================================================ #

function factorisation()
    fig, ax = canvas(820, 300)

    c = 16                    # side of a cell
    g = 4                     # gap between two cells
    w = 3*c + 2*g             # width of a matrix with three columns
    h = 150                   # height of a matrix of three states
    y = 158                   # height of the middle of the matrices

    # a rectangle with its bottom left corner at (x, y)
    box!(x, y, w, h; color, strokecolor=:transparent, linestyle=:solid) =
        poly!(ax, Rect2f(x, y, w, h); color, strokecolor, strokewidth=1.2, linestyle)

    # a matrix of three by three cells with its left side at x, each cell coloured by `fill`,
    # which takes the row and the column and returns a colour or nothing
    function square!(x, fill)
        for i in 1:3, j in 1:3
            col = fill(i, j)
            yᵢ  = y + w/2 - i*c - (i - 1)*g
            xⱼ  = x + (j - 1)*(c + g)
            if col === nothing
                box!(xⱼ, yᵢ, c, c; color=:white, strokecolor=GRID)
            else
                box!(xⱼ, yᵢ, c, c; color=col)
            end

            # the entries equal to one, on the diagonal of L and of its inverse
            col == (INK2, 0.9) &&
                label!(ax, Point2f(xⱼ + c/2, yᵢ + c/2), "1"; fontsize=11, color=:white, font=:bold)
        end
    end

    # a matrix of three states with its left side at x, the column j coloured by `fill`
    function tall!(x, fill; linestyle=:solid)
        for j in 1:3
            col = fill(j)
            box!(x + (j - 1)*(c + g), y - h/2, c, h; color=col[1], strokecolor=col[2], linestyle)
        end
    end

    gray  = (MUTED, 0.45)
    upper = (i, j) -> i <= j ? gray : nothing
    lower = (i, j) -> i > j ? gray : i == j ? (INK2, 0.9) : nothing
    bycol = (i, j) -> i <= j ? SLOTS[j] : nothing
    unit  = (i, j) -> i < j ? SLOTS[j] : i == j ? (INK2, 0.9) : nothing

    # ---- the inner products and their factorisation ----
    x = 34
    square!(x, (i, j) -> gray)
    label!(ax, Point2f(x + w + 16, y), L"="; fontsize=20)
    square!(x + w + 32, lower)
    square!(x + 2w + 44, upper)
    label!(ax, Point2f(x + w/2, y + w/2 + 18), L"P")
    label!(ax, Point2f(x + 1.5w + 32, y + w/2 + 18), L"L")
    label!(ax, Point2f(x + 2.5w + 44, y + w/2 + 18), L"U")
    label!(ax, Point2f(x + 1.5w + 22, y - w/2 - 30), L"P_{ij} = \langle \psi_i^+, \psi_j^- \rangle";
           fontsize=16)
    label!(ax, Point2f(x + 1.5w + 22, 30), "inner products of the two sets"; fontsize=13, color=INK2)

    # ---- the CLVs ----
    x = 318
    tall!(x, j -> (SLOTS[j], :transparent))
    label!(ax, Point2f(x + w + 16, y), L"="; fontsize=20)
    tall!(x + w + 32, j -> (gray, :transparent))
    square!(x + 2w + 44, bycol)
    label!(ax, Point2f(x + w/2, y + h/2 + 18), L"\tilde{\Gamma}")
    label!(ax, Point2f(x + w/2, y - h/2 - 14), L"N"; fontsize=15, color=INK2)
    label!(ax, Point2f(x - 14, y), L"n"; fontsize=15, color=INK2)
    label!(ax, Point2f(x + 1.5w + 32, y + h/2 + 18), L"\Psi^-")
    label!(ax, Point2f(x + 2.5w + 44, y + w/2 + 18), L"U^{-1}")
    label!(ax, Point2f(x + 1.5w + 22, 30), "CLVs, from the tangent set"; fontsize=13, color=INK2)

    # ---- the adjoint CLVs ----
    x = 590
    tall!(x, j -> ((SLOTS[j], 0.35), SLOTS[j]); linestyle=DASH)
    label!(ax, Point2f(x + w + 16, y), L"="; fontsize=20)
    tall!(x + w + 32, j -> ((MUTED, 0.15), MUTED); linestyle=DASH)
    square!(x + 2w + 44, unit)
    label!(ax, Point2f(x + w/2, y + h/2 + 18), L"\tilde{\Theta}")
    label!(ax, Point2f(x + 1.5w + 32, y + h/2 + 18), L"\Psi^+")
    label!(ax, Point2f(x + 2.5w + 44, y + w/2 + 18), L"L^{-\top}")
    label!(ax, Point2f(x + 1.5w + 22, 30), "adjoint CLVs, from the adjoint set"; fontsize=13, color=INK2)

    save(joinpath(ASSETS, "factorisation.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Algorithm: the two sweeps along the trajectory                               #
# ============================================================================ #

function algorithm()
    fig, ax = canvas(820, 410)

    # ---- times, in units of the figure ----
    x₀ = 176.0                # time zero
    x₁ = 344.0                # start of the window, t_B
    x₂ = 614.0                # end of the window, t_C
    x₃ = 782.0                # end of the second spin-up
    M  = 3                    # chunks of a spin-up
    K  = 6                    # sampling intervals of the window
    Δ  = (x₁ - x₀)/M
    xs = range(x₁, x₂; length=K + 1)

    # ---- heights of the rows ----
    y₁ = 322                  # trajectory, forward sweep
    y₂ = 272                  # tangent set
    y₃ = 202                  # trajectory, backward sweep
    y₄ = 152                  # adjoint set
    y₅ = 86                   # the two families
    y₆ = 40                   # axis of time

    # the number of a step of the algorithm, in a disc
    function step!(p, k)
        scatter!(ax, [p]; color=INK2, markersize=19)
        label!(ax, p, string(k); fontsize=12, color=:white, font=:bold)
    end

    # ---- the window and the two spin-ups ----
    poly!(ax, Rect2f(x₁, y₅ - 26, x₂ - x₁, 388 - y₅ + 26); color=WASH, strokewidth=0)
    for (x, s) in (((x₀ + x₁)/2, "spin-up"), ((x₁ + x₂)/2, "window"), ((x₂ + x₃)/2, "spin-up"))
        label!(ax, Point2f(x, 374), s; fontsize=13, color=INK2)
    end

    # ---- the samples of the two sets are combined at each time ----
    for x in xs
        segment!(ax, Point2f(x, y₂), Point2f(x, y₅); color=AXIS)
    end

    # ---- forward sweep: the trajectory, in chunks, and the tangent set ----
    for k in 1:M
        arrow!(ax, Point2f(x₀ + (k - 1)*Δ, y₁), Point2f(x₀ + k*Δ - 2, y₁); color=MUTED, linewidth=2, head=9)
    end
    arrow!(ax, Point2f(x₁, y₁), Point2f(x₂ - 2, y₁); color=MUTED, linewidth=2, head=9)
    for k in 1:M-1
        arrow!(ax, Point2f(x₂ + (k - 1)*Δ, y₁), Point2f(x₂ + k*Δ - 6, y₁); color=MUTED, linewidth=2, head=9)
    end
    scatter!(ax, Point2f[(x₂ + k*Δ, y₁) for k in 0:M-1]; color=:white, markersize=10,
             strokecolor=INK2, strokewidth=1.6)
    arrow!(ax, Point2f(x₀, y₂), Point2f(x₂ + 24, y₂); color=INK)
    scatter!(ax, Point2f[(x, y₂) for x in xs]; color=INK, markersize=12, strokecolor=:white,
             strokewidth=1.5)

    step!(Point2f((x₀ + x₁)/2, y₁ + 20), 1)
    step!(Point2f((x₁ + x₂)/2, y₁ + 20), 2)
    step!(Point2f(x₂ + Δ, y₁ + 20), 3)

    # ---- backward sweep: the trajectory again, from the checkpoints, and the adjoint set ----
    arrow!(ax, Point2f(x₁, y₃), Point2f(x₂ - 6, y₃); color=MUTED, linewidth=2, head=9)
    for k in 1:M
        arrow!(ax, Point2f(x₂ + (k - 1)*Δ, y₃), Point2f(x₂ + k*Δ - 6, y₃); color=MUTED, linewidth=2, head=9)
    end
    scatter!(ax, Point2f[(x₂ + k*Δ, y₃) for k in 0:M-1]; color=:white, markersize=10,
             strokecolor=INK2, strokewidth=1.6)
    arrow!(ax, Point2f(x₃, y₄), Point2f(x₁ - 24, y₄); color=INK, linestyle=DASH)
    scatter!(ax, Point2f[(x, y₄) for x in xs]; color=:white, markersize=10, strokecolor=INK, strokewidth=2)

    step!(Point2f((x₂ + x₃)/2, y₃ + 20), 4)
    step!(Point2f((x₁ + x₂)/2, y₃ + 20), 5)

    # ---- the two families ----
    for x in xs
        scatter!(ax, Point2f[(x, y₅)]; color=INK, marker=:rect, markersize=11, strokecolor=:white,
                 strokewidth=2)
    end
    step!(Point2f(x₂ + 26, y₅), 6)
    label!(ax, Point2f(x₂ + 44, y₅), L"P = LU"; align=(:left, :center), fontsize=16)

    # ---- names of the rows ----
    for (y, s) in ((y₁, "trajectory"), (y₂, "tangent set"), (y₃, "trajectory, again"),
                   (y₄, "adjoint set"), (y₅, "CLVs and adjoint CLVs"))
        label!(ax, Point2f(16, y), s; align=(:left, :center), fontsize=13, color=INK2)
    end
    label!(ax, Point2f(x₀ + 6, y₂ + 15), "random"; align=(:left, :center), fontsize=12, color=MUTED)
    label!(ax, Point2f(x₃ - 6, y₄ + 15), "random"; align=(:right, :center), fontsize=12, color=MUTED)
    label!(ax, Point2f(x₀ + 64, y₂ - 16), L"\psi^-"; fontsize=16)
    label!(ax, Point2f(x₃ - 64, y₄ - 16), L"\psi^+"; fontsize=16)

    # ---- the axis of time ----
    arrow!(ax, Point2f(x₀ - 14, y₆), Point2f(x₃ + 26, y₆); color=AXIS, linewidth=1.5, head=9)
    for (x, s) in ((x₀, L"0"), (x₁, L"t_B"), (x₂, L"t_C"), (x₃, L"t_C + t_\mathrm{spin}"))
        segment!(ax, Point2f(x, y₆ - 5), Point2f(x, y₆ + 5); color=MUTED)
        label!(ax, Point2f(x, y₆ - 20), s; fontsize=15)
    end
    label!(ax, Point2f(16, y₆), "time"; align=(:left, :center), fontsize=13, color=INK2)

    save(joinpath(ASSETS, "algorithm.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Sampling: the samples of the two sets must fall at the same times            #
# ============================================================================ #

function sampling()
    fig, ax = canvas(820, 270)

    x₁ = 196.0                # start of the window, t_B
    h  = 132.0                # sampling interval, of four time steps
    δ  = 17                   # distance of the samples from the axis of time

    for (y, K, name, verdict) in ((200, 3.0, "(a)", "same times"), (72, 3.5, "(b)", "ArgumentError"))
        x₂ = x₁ + K*h

        # ---- the window and its time steps ----
        segment!(ax, Point2f(x₁, y), Point2f(x₂, y); color=AXIS, linewidth=1.5)
        for x in x₁:h/4:x₂
            segment!(ax, Point2f(x, y - 3), Point2f(x, y + 3); color=AXIS)
        end
        label!(ax, Point2f(x₁ - 16, y), L"t_B"; align=(:right, :center), fontsize=15)
        label!(ax, Point2f(x₂ + 16, y), L"t_C"; align=(:left, :center), fontsize=15)

        # ---- samples of the tangent set, counted from the start, and the last state ----
        tangent = unique([collect(x₁:h:x₂); x₂])

        # ---- samples of the adjoint set, counted from the end, and the last state ----
        adjoint = sort(unique([collect(x₂:-h:x₁); x₁]))

        # ---- the samples are combined in pairs, in the order of time ----
        for (a, b) in zip(tangent, adjoint)
            segment!(ax, Point2f(a, y + δ), Point2f(b, y - δ); color = a == b ? MUTED : "#d03b3b",
                     linewidth = a == b ? 1 : 1.6)
        end
        scatter!(ax, Point2f[(x, y + δ) for x in tangent]; color=INK, markersize=12,
                 strokecolor=:white, strokewidth=1.5)
        scatter!(ax, Point2f[(x, y - δ) for x in adjoint]; color=:white, markersize=10,
                 strokecolor=INK, strokewidth=2)

        # ---- the direction in which the time steps are counted ----
        arrow!(ax, Point2f(x₁, y + δ + 18), Point2f(x₁ + 56, y + δ + 18); color=INK2, linewidth=1.3, head=8)
        arrow!(ax, Point2f(x₂, y - δ - 18), Point2f(x₂ - 56, y - δ - 18); color=INK2, linewidth=1.3,
               head=8, linestyle=DASH)

        label!(ax, Point2f(x₁ - 16, y + δ + 14), "tangent set"; align=(:right, :center), fontsize=13,
               color=INK2)
        label!(ax, Point2f(x₁ - 16, y - δ - 14), "adjoint set"; align=(:right, :center), fontsize=13,
               color=INK2)
        label!(ax, Point2f(700, y + δ + 18), verdict; align=(:left, :center), fontsize=13, color=INK2)
        panellabel!(ax, Point2f(16, y + 56), name)
    end

    save(joinpath(ASSETS, "sampling.png"), fig; px_per_unit=3)
    return nothing
end


# ============================================================================ #
# Draw the diagrams named on the command line, or all of them                  #
# ============================================================================ #

for diagram in (covariance, duality, transversality, directions,
                subspaces, factorisation, algorithm, sampling)
    isempty(ARGS) || string(diagram) in ARGS || continue
    diagram()
    println("written assets/$(diagram).png")
end
