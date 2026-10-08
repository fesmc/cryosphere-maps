# Antarctica A0 landscape poster map.
# Usage: julia --project=. scripts/antarctica.jl [velocity|surface|bed] [dark] [nocontours] [cmap=<name>] [draft]
# Requires the prepared grid from steps 0-1 (see README).

include("common.jl")

const CREDITS_ANT = "Data: bed and surface topography from BedMachine Antarctica v4 (Morlighem et al., 2020); " *
    "surface ice velocity from MEaSUREs InSAR v2 (Rignot et al., 2011); drainage basins from IMBIE 2 " *
    "(Rignot & Mouginot); place names from the SCAR Composite Gazetteer of Antarctica. " *
    "Polar stereographic projection (71°S, 0°E), $(round(Int, 1000GRIDS["antarctica"].dx)) m grid."

load_antarctica(; stride=1) = load_prepared("antarctica"; zfac_ice=40, stride)

function plot_antarctica(d, style; labels=nothing, maxtier=1, scale=1.85, ocean=:light, contours=true,
                         velcmap=CS_VEL)
    xmap = (-3040.0, 3040.0); ymap = (-2600.0, 2500.0)
    pad  = 50; panel = 520.0
    H    = A0_LANDSCAPE[2] - 2pad
    kmpp = (ymap[2] - ymap[1])/H
    W    = (xmap[2] - xmap[1])/kmpp

    r   = d.r
    img = compose(style, r; srflim=(0, 4100), bedlim=(-2500, 3000), ocean, velcmap)

    fig = Figure(size=A0_LANDSCAPE, figure_padding=pad, backgroundcolor=:white, fontsize=24)
    ax  = Axis(fig[1, 1]; width=W, height=H, limits=(xmap, ymap), backgroundcolor=:white)
    hidedecorations!(ax); hidespines!(ax)

    ix = findall(xi -> xmap[1] <= xi <= xmap[2], r.x); iy = findall(yi -> ymap[1] <= yi <= ymap[2], r.y)
    h = r.dx/2
    image!(ax, (r.x[ix[1]] - h, r.x[ix[end]] + h), (r.y[iy[1]] - h, r.y[iy[end]] + h), img[ix, iy]; interpolate=true)
    xs = r.x[ix]; ys = r.y[iy]

    if contours || style == :surface
        elevation_contours!(ax, xs, ys, r.zs[ix, iy], d.grounded[ix, iy], 500:500:4000)
    end
    # ice front and grounding line
    maskoutline!(ax, xs, ys, r.ice[ix, iy]; σ=1.5, color=(:black, 0.75), linewidth=1.5)
    maskoutline!(ax, xs, ys, d.grounded[ix, iy]; σ=1.5, mask=r.ice[ix, iy], color=(:gray25, 0.8), linewidth=1.2)
    # drainage divides over grounded ice, on the coarser basin grid
    basin_divides!(ax, d.xb, d.yb, d.basin, d.grounded_b; σ=1.5, color=(:black, 0.6), lw=2.5)

    graticule!(ax, PROJ_ANT, -80:10:-60, -180:30:150; latrange=(-88, -55), color=(:gray20, 0.4), lw=1.0)
    scalebar!(ax, xmap[1] + 150, ymap[1] + 150, 1000.0; fontsize=26)

    if labels !== nothing
        draw_labels!(ax, labels, PROJ_ANT; layout=:coastal, kmpp, maxtier, scale,
                     icemask=r.ice, shelfmask=r.shelf, groundedmask=d.grounded, x=r.x, y=r.y, offset=120.0,
                     limits=(xmap..., ymap...), seacolor=seacolor(style, ocean))
    end

    side = fig[1, 2] = GridLayout(width=panel)
    top  = side[1, 1] = GridLayout(valign=:top, tellheight=false)
    Label(top[1, 1], "Antarctic\nIce Sheet"; fontsize=96, font=:bold, halign=:left, justification=:left, lineheight=0.95)
    Label(top[2, 1], "Surface ice velocity, drainage\ndivides and place names"; fontsize=36, color=:gray30,
          halign=:left, justification=:left)
    style_colorbar!(top[3, 1], style; velcmap, vertical=false, width=panel - 40, height=32, labelsize=32, ticklabelsize=28,
                    halign=:left)
    symbol_legend!(top[4, 1]; scale, labelsize=28, halign=:left,
                   extra=[(LineElement(color=(:gray25, 0.8), linewidth=1.2), "Grounding line")])
    rowgap!(top, 1, 20); rowgap!(top, 2, 100); rowgap!(top, 3, 80)
    Label(side[2, 1], CREDITS_ANT; fontsize=22, color=:gray30, word_wrap=true, width=panel, justification=:left,
          halign=:left, valign=:bottom, tellheight=false)
    colgap!(fig.layout, 60)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    o = parse_args(ARGS)
    d = load_antarctica(; o.stride)
    labs = read_labels(joinpath(ROOT, "data", "labels_antarctica.csv"))
    fig = plot_antarctica(d, o.style; labels=labs, maxtier=1, o.ocean, o.contours, o.velcmap)
    save_poster(fig, joinpath(ROOT, "plots", "antarctica_A0_$(o.style)$(o.tag)"); o.draft)
end
