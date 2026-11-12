# Colours and theme of the figures of the README, shared by the scripts of this folder.
#
#     include(joinpath(@__DIR__, "style.jl"))

using CairoMakie

# The first three slots of a categorical palette validated for colour vision deficiency (all
# pairs). The first, second and third vector of a family keep their colour in every figure. On
# the Lorenz equations these are the unstable, neutral and stable directions.
const UNSTABLE = "#eb6834"
const NEUTRAL  = "#1baf7a"
const STABLE   = "#2a78d6"
const SLOTS    = (UNSTABLE, NEUTRAL, STABLE)
const INK      = "#0b0b0b"
const INK2     = "#52514e"
const MUTED    = "#898781"
const GRID     = "#e1e0d9"
const AXIS     = "#c3c2b7"
const WASH     = "#f4f3ee"

# A sequential map for magnitudes, one hue from light to dark.
const SEQUENTIAL = cgrad(["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])

set_theme!(Theme(fontsize        = 14,
                 backgroundcolor = :white,
                 figure_padding  = (8, 16, 8, 8),
                 Axis   = (backgroundcolor    = :white,
                           xgridcolor         = GRID,
                           ygridcolor         = GRID,
                           xminorgridvisible  = false,
                           yminorgridvisible  = false,
                           leftspinecolor     = AXIS,
                           bottomspinecolor   = AXIS,
                           topspinevisible    = false,
                           rightspinevisible  = false,
                           xtickcolor         = AXIS,
                           ytickcolor         = AXIS,
                           xticklabelcolor    = INK2,
                           yticklabelcolor    = INK2,
                           xlabelcolor        = INK,
                           ylabelcolor        = INK,
                           titlecolor         = INK),
                 Legend = (framevisible = false,
                           labelcolor   = INK2,
                           patchsize    = (22, 12)),
                 Lines  = (linewidth = 2,),
                 Scatter = (markersize = 9,)))

# the label of a panel, in its top left corner
panel!(ax, s) = text!(ax, 0.02, 0.98; text=s, space=:relative, align=(:left, :top),
                      font=:bold, color=INK, fontsize=15)
