# Draws assets/logo.svg: a portrait of Aleksandr Lyapunov wearing pixelated sunglasses, in front
# of two trajectories of the Lorenz equations that start close to each other, over a landscape
# that recedes to the horizon, beside the CLVs.jl wordmark.
#
#     julia assets/make_logo.jl
#
# The portrait, assets/lyapunov.png, is a photograph in the public domain, with the background
# removed. It was published in 1908 by I. P. Osipov and D. I. Bagalei in the history of the
# Faculty of Physics and Mathematics of the University of Kharkiv, and is taken from Wikimedia
# Commons, File:Aleksandr Lyapunov.jpg.

using Base64
using Printf


# ---- the frame of the portrait, in the units of the svg ----

X₀ = 58.0    # left side
Y₀ = 48.0    # top side
W  = 360.0   # width
H  = 204.0   # height

# colours of Julia, and the colour of the background of the frame
PURPLE = (149, 88, 178)
BLUE   = (64, 99, 216)
PAPER  = (251, 250, 248)

# the colour between a and b, at the fraction t from a
mix(a, b, t) = @sprintf("#%02x%02x%02x", (round(Int, (1 - t)*a[i] + t*b[i]) for i in 1:3)...)

# the path through a list of points
path(points) = "M" * join((@sprintf("%.1f %.1f", p...) for p in points), " L")


# ---- the landscape ----

# The eye is at the height EYE over the floor of a valley that runs to the horizon. A point at x
# across, at the depth y from the eye and at the height z is drawn at (X, Y).
XC    = X₀ + W/2   # vanishing point, behind the portrait
YH    = Y₀ + 138   # horizon
FOCAL = 180.0      # focal length
EYE   = 2.2        # height of the eye

project(x, y, z) = (XC + FOCAL*x/y, YH - FOCAL*(z - EYE)/y)

# height of the ground: a valley between low hills, with a range of mountains far away
smooth(s) = s <= 0 ? 0.0 : s >= 1 ? 1.0 : s^2*(3 - 2s)

function ground(x, y)
    hills  = 1.4*(1 - exp(-(x/4.5)^2))*(0.6 + 0.4*sin(0.55y + 0.12x))
    range  = 26*smooth((y - 18)/36)*(0.15 + 0.85*(1 - exp(-(x/24)^2)))*(0.62 + 0.38*sin(0.19x + 1.7)*cos(0.11y))
    ripple = 0.3*sin(1.3x + 0.7y)*sin(0.9y - 0.4x)
    return hills + range + ripple
end

# Lines across the valley at depths in geometric progression, and lines along it every 1.5. Each
# strip between two lines across is filled before the nearer ones, which hide what lies behind.
depths    = [3.9*1.085^j for j in 0:34]
landscape = String[]
for j in length(depths)-1:-1:1
    near = depths[j]
    far  = depths[j+1]
    t    = (j - 1)/(length(depths) - 2)   # 0 for the nearest strip, 1 for the farthest
    half = (W/2 + 40)*far/FOCAL           # half width of the strip that is in view

    # the strip, with the nearer line from left to right and the farther one back
    xs    = range(-half, half; length=81)
    front = [project(x, near, ground(x, near)) for x in xs]
    back  = [project(x, far, ground(x, far)) for x in reverse(xs)]
    push!(landscape, @sprintf("<path d=\"%s Z\" fill=\"%s\"/>",
                              path(vcat(front, back)), mix((236, 230, 247), PAPER, t^0.6)))

    # the lines of the grid on the strip, fainter with the distance
    along = join((path([project(x, far, ground(x, far)), project(x, near, ground(x, near))])
                  for x in -1.5*floor(half/1.5):1.5:half), " ")
    push!(landscape, @sprintf("<path d=\"%s %s\" fill=\"none\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-opacity=\"%.2f\"/>",
                              path(front), along, mix(BLUE, PURPLE, t), 0.9 - 0.4*t, 0.08 + 0.5*(1 - t)^1.4))
end


# ---- two trajectories of the Lorenz equations ----

lorenz(u) = [10*(u[2] - u[1]), 28*u[1] - u[2] - u[1]*u[3], -8/3*u[3] + u[1]*u[2]]

# one step of size h of the classical Runge-Kutta scheme
function rk4(u, h)
    k₁ = lorenz(u)
    k₂ = lorenz(u .+ h/2 .* k₁)
    k₃ = lorenz(u .+ h/2 .* k₂)
    k₄ = lorenz(u .+ h .* k₃)
    return u .+ h/6 .* (k₁ .+ 2 .* k₂ .+ 2 .* k₃ .+ k₄)
end

# The projection on the plane (x, z) of n steps from u, scaled to the frame, split by the
# distance from the viewer, who looks along y, into LAYERS paths from the nearest to the farthest.
LAYERS = 6

function trajectory(u, n, h)
    layers = [String[] for l in 1:LAYERS]
    last   = 0
    for i in 1:n
        v = rk4(u, h)
        l = clamp(1 + floor(Int, LAYERS*((u[2] + v[2])/2 + 25)/50), 1, LAYERS)
        p = @sprintf("%.1f %.1f", X₀ + W/2 + 8.4*u[1], Y₀ + H - 4 - 4.1*(u[3] - 1))
        q = @sprintf("%.1f %.1f", X₀ + W/2 + 8.4*v[1], Y₀ + H - 4 - 4.1*(v[3] - 1))
        push!(layers[l], l == last ? "L$q" : "M$p L$q")
        last = l
        u    = v
    end
    return join.(layers, " ")
end

# a point on the attractor, and one at a distance of 0.01 from it
u = [1.0, 1.0, 20.0]
for i in 1:4000
    global u = rk4(u, 0.005)
end
first  = trajectory(u, 2200, 0.005)
second = trajectory(u .+ [0.01, 0.0, 0.0], 2200, 0.005)

# the farthest layers first, thinner and fainter than the nearest
attractor = String[]
for l in LAYERS:-1:1, (d, colour) in ((first[l], "#9558b2"), (second[l], "#4063d8"))
    t = (l - 1)/(LAYERS - 1)
    push!(attractor, @sprintf("<path d=\"%s\" fill=\"none\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-opacity=\"%.2f\" stroke-linejoin=\"round\"/>",
                              d, colour, 1.9 - 1.25*t, 0.9 - 0.5*t))
end


# ---- the portrait ----

# The photograph, of 600 by 783 pixels, is scaled by S and its top left corner is at (PX, PY).
# The face is in the middle of the frame. The bust fades at its sides, where the photograph
# ends, and a mist covers the bottom of the picture.
S  = 0.32
PX = X₀ + W/2 - S*303
PY = Y₀ + 12 - S*90

photo = base64encode(read(joinpath(@__DIR__, "lyapunov.png")))


# ---- the sunglasses ----

# The pixels, X for black and o for white. The two lenses are 15 pixels apart, and are drawn
# slightly larger than the distance between the eyes would give.
PIXELS = ["XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
          "  XXXXXXXXXXXXXXXXXXXXXXXXXXXXX  ",
          "   XXoXoXXXXXXX   XXoXoXXXXXXX   ",
          "   XXXoXoXXXXXX   XXXoXoXXXXXX   ",
          "    XXXoXoXXXX     XXXoXoXXXX    ",
          "     XXXXXXXX       XXXXXXXX     ",
          "      XXXXXX         XXXXXX      "]

# The eyes in the photograph, in pixels, set the position, the size and the tilt of the glasses.
eyeₗ = (257.5, 275.0)
eyeᵣ = (347.5, 269.0)
unit = 1.1 * hypot(eyeᵣ[1] - eyeₗ[1], eyeᵣ[2] - eyeₗ[2]) / 15
tilt = atand(eyeᵣ[2] - eyeₗ[2], eyeᵣ[1] - eyeₗ[1])

pixels = String[]
for (i, row) in enumerate(PIXELS), (j, c) in enumerate(row)
    c == ' ' && continue
    fill = c == 'X' ? "#111318" : "#ffffff"
    push!(pixels, @sprintf("<rect x=\"%d\" y=\"%d\" width=\"1.04\" height=\"1.04\" fill=\"%s\"/>",
                           j - 1, i - 1, fill))
end

# the middle of the glasses, between the lenses and at the height of the eyes, goes between the eyes
glasses = @sprintf("<g transform=\"translate(%.2f %.2f) scale(%.3f) translate(%.2f %.2f) rotate(%.2f) scale(%.3f) translate(%.1f %.1f)\">",
                   PX, PY, S,
                   (eyeₗ[1] + eyeᵣ[1])/2, (eyeₗ[2] + eyeᵣ[2])/2,
                   tilt, unit,
                   -length(PIXELS[1])/2, -3.6) * join(pixels) * "</g>"


# ---- the svg ----

svg = """
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 980 300" role="img" aria-labelledby="title desc">
  <title id="title">CLVs.jl logo</title>
  <desc id="desc">A portrait of Aleksandr Lyapunov wearing pixelated sunglasses, in front of two trajectories of the Lorenz equations in Julia colors over a landscape that recedes to the horizon, and the CLVs.jl wordmark.</desc>
  <defs>
    <linearGradient id="panel" x1="0" x2="1" y1="0" y2="1">
      <stop offset="0" stop-color="#ffffff"/>
      <stop offset="1" stop-color="#f6f8fb"/>
    </linearGradient>
    <filter id="shadow" x="-10%" y="-20%" width="120%" height="145%">
      <feDropShadow dx="0" dy="8" stdDeviation="10" flood-color="#172033" flood-opacity="0.14"/>
    </filter>
    <clipPath id="portrait-clip">
      <rect x="58" y="48" width="360" height="204" rx="26"/>
    </clipPath>
    <linearGradient id="sky" x1="0" x2="0" y1="0" y2="1">
      <stop offset="0" stop-color="#e3dcf4"/>
      <stop offset="0.68" stop-color="#fbf7f4"/>
    </linearGradient>
    <radialGradient id="glow" gradientUnits="userSpaceOnUse" cx="$(XC)" cy="$(YH - 4)" r="170">
      <stop offset="0" stop-color="#f5b9a0" stop-opacity="0.75"/>
      <stop offset="0.45" stop-color="#f3c9d6" stop-opacity="0.35"/>
      <stop offset="1" stop-color="#f3c9d6" stop-opacity="0"/>
    </radialGradient>
    <linearGradient id="mist" x1="0" x2="0" y1="0" y2="1">
      <stop offset="0.76" stop-color="#fbfaf8" stop-opacity="0"/>
      <stop offset="1" stop-color="#fbfaf8" stop-opacity="0.92"/>
    </linearGradient>
    <linearGradient id="fade-left" x1="0" x2="1" y1="0" y2="0">
      <stop offset="0.2" stop-color="#000000"/>
      <stop offset="1" stop-color="#000000" stop-opacity="0"/>
    </linearGradient>
    <linearGradient id="fade-right" x1="1" x2="0" y1="0" y2="0">
      <stop offset="0.2" stop-color="#000000"/>
      <stop offset="1" stop-color="#000000" stop-opacity="0"/>
    </linearGradient>
    <mask id="bust" maskUnits="userSpaceOnUse" x="58" y="48" width="360" height="204">
      <rect x="58" y="48" width="360" height="204" fill="#ffffff"/>
      <rect x="$(@sprintf("%.1f", PX))" y="48" width="30" height="204" fill="url(#fade-left)"/>
      <rect x="$(@sprintf("%.1f", PX + 600S - 30))" y="48" width="31" height="204" fill="url(#fade-right)"/>
    </mask>
  </defs>

  <rect x="24" y="24" width="932" height="252" rx="28" fill="url(#panel)" filter="url(#shadow)"/>

  <g clip-path="url(#portrait-clip)">
    <rect x="58" y="48" width="360" height="204" fill="url(#sky)"/>
    <rect x="58" y="48" width="360" height="204" fill="url(#glow)"/>
    $(join(landscape, "\n    "))
    $(join(attractor, "\n    "))
    <g mask="url(#bust)">
      <image x="$(@sprintf("%.2f", PX))" y="$(@sprintf("%.2f", PY))" width="$(@sprintf("%.2f", 600S))" height="$(@sprintf("%.2f", 783S))" href="data:image/png;base64,$photo"/>
    </g>
    $glasses
    <rect x="58" y="48" width="360" height="204" fill="url(#mist)"/>
  </g>
  <rect x="58" y="48" width="360" height="204" rx="26" fill="none" stroke="#d7dce6" stroke-width="2"/>

  <g transform="translate(462 70)" font-family="Inter, Avenir, Helvetica, Arial, sans-serif">
    <text x="0" y="62" font-size="64" font-weight="800" fill="#242936">CLVs.jl</text>
    <text x="2" y="104" font-size="20" font-weight="700" fill="#5c6475">covariant Lyapunov vectors and their adjoints</text>
    <text x="2" y="133" font-size="20" font-weight="600" fill="#5c6475">for trajectories of large dynamical systems</text>
  </g>
</svg>
"""

write(joinpath(@__DIR__, "logo.svg"), svg)
println("written assets/logo.svg")
