# Greenland A0 portrait poster map.
# Usage: julia --project=. scripts/greenland.jl [velocity|surface|bed ...]
# Requires the prepared grid from steps 0-1 (see README).

include("common.jl")

const CREDITS_GRL = "Data: bed and surface topography from BedMachine Greenland v6 (Morlighem et al., 2017); " *
    "surface ice velocity from the MEaSUREs multi-year mosaic (Joughin et al., 2018); drainage regions from " *
    "IMBIE 2 (Rignot & Mouginot); outside the BedMachine domain, topography from ETOPO 2022 (NOAA NCEI, 2022) " *
    "and glacier outlines from RGI 6.0 (RGI Consortium, 2017); place names from GeoNames. " *
    "Polar stereographic projection (70°N, 45°W), $(round(Int, 1000GRIDS["greenland"].dx)) m grid."

load_greenland() = load_prepared("greenland"; zfac_ice=40)

function plot_greenland(d, style; labels=nothing, maxtier=1, scale=1.85)
    # Map window (km) fills the page width; beyond the BedMachine domain the
    # topography/bathymetry is ETOPO 2022 and glacier ice is from RGI 6.0.
    pad  = 70
    W    = A0_PORTRAIT[1] - 2pad
    H    = 2790.0
    ymap = (-3420.0, -600.0)
    kmpp = (ymap[2] - ymap[1])/H
    xc   = 120.0                                  # centre of Greenland
    xl   = (xc - W*kmpp/2, xc + W*kmpp/2)

    r   = d.r
    img = compose(style, r; srflim=(0, 3300), bedlim=(-1500, 3000))

    fig = Figure(size=A0_PORTRAIT, figure_padding=pad, backgroundcolor=:white, fontsize=24)
    head = fig[1, 1] = GridLayout()
    Label(head[1, 1], "Greenland Ice Sheet"; fontsize=120, font=:bold)
    Label(head[2, 1], "Surface ice velocity, drainage regions and place names"; fontsize=44, color=:gray30)
    rowgap!(head, 10)

    ax = Axis(fig[2, 1]; width=W, height=H, limits=(xl, ymap), backgroundcolor=:white)
    hidedecorations!(ax); hidespines!(ax)

    ix = findall(xi -> xl[1] <= xi <= xl[2], r.x); iy = findall(yi -> ymap[1] <= yi <= ymap[2], r.y)
    h = r.dx/2
    image!(ax, (r.x[ix[1]] - h, r.x[ix[end]] + h), (r.y[iy[1]] - h, r.y[iy[end]] + h), img[ix, iy]; interpolate=true)
    lines!(ax, Rect(xl[1], ymap[1], xl[2] - xl[1], ymap[2] - ymap[1]); color=:black, linewidth=1.5)

    xs = r.x[ix]; ys = r.y[iy]
    if style == :surface
        contourlines!(ax, xs, ys, r.zs[ix, iy], 500:500:3000; mask=r.ice[ix, iy],
                      color=(:steelblue4, 0.45), linewidth=0.9)
    end
    maskoutline!(ax, xs, ys, r.ice[ix, iy]; σ=1.5, color=(:black, 0.75), linewidth=1.5)
    basin_divides!(ax, d.xb, d.yb, d.basin, d.ice_b; σ=1.2, color=(:black, 0.6), lw=2.5)
    regions = basin_labels!(ax, d.xb, d.yb, d.basin, d.basin_names, d.ice_b; kmpp, fontsize=44, color=(:gray20, 0.8))

    graticule!(ax, PROJ_GRL, 60:5:80, -80:10:0; latrange=(58, 84), color=(:gray20, 0.4), lw=1.0,
               inside=p -> xl[1] <= p[1] <= xl[2] && ymap[1] <= p[2] <= ymap[2])
    scalebar!(ax, 700.0, ymap[1] + 60, 400.0; fontsize=26, color=:white)

    if labels !== nothing
        draw_labels!(ax, labels, PROJ_GRL; layout=:coastal, kmpp, maxtier, scale, obstacles=regions,
                     icemask=r.ice, x=r.x, y=r.y, offset=50.0, maxlead=300.0, tmax=120.0,
                     seacolor=(style == :bed ? :gray20 : :white), limits=(xl..., ymap...))
    end

    foot = fig[3, 1] = GridLayout()
    style_colorbar!(foot[1, 1], style; vertical=false, width=600, height=32, labelsize=32, ticklabelsize=28)
    symbol_legend!(foot[1, 2]; scale, domes=false, labelsize=28, nbanks=2, tellheight=true)
    colgap!(foot, 150)
    Label(fig[4, 1], CREDITS_GRL; fontsize=22, color=:gray30, word_wrap=true, width=W - 200, justification=:left)
    rowgap!(fig.layout, 30)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    styles = isempty(ARGS) ? [:velocity] : Symbol.(ARGS)
    d = load_greenland()
    labs = read_labels(joinpath(ROOT, "data", "labels_greenland.csv"))
    for s in styles
        save_poster(plot_greenland(d, s; labels=labs, maxtier=1), joinpath(ROOT, "plots", "greenland_A0_$(s)"))
    end
end
