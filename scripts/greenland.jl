# Greenland A0 portrait poster map.
# Usage: julia --project=. scripts/greenland.jl [velocity|surface|bed] [dark] [nocontours] [cmap=<name>]
#            [tier=2] [add=<name>;...] [list[=<text>]] [full]
# Options: see parse_args in common.jl and the README.

include("common.jl")

const CREDITS_GRL = "Data: bed and surface topography from BedMachine Greenland v6 (Morlighem et al., 2017); " *
    "surface ice velocity from the MEaSUREs multi-year mosaic (Joughin et al., 2018); drainage regions from " *
    "IMBIE 2 (Rignot & Mouginot); outside the BedMachine domain, topography from ETOPO 2022 (NOAA NCEI, 2022) " *
    "and glacier outlines from RGI 6.0 (RGI Consortium, 2017); median sea-ice edge 1981–2010 from the Sea Ice Index " *
    "v4 (Fetterer et al., 2017); place names from GeoNames. Key numbers for the ice sheet within the IMBIE 2 regions, " *
    "computed on the map grid; sea-level equivalent of the ice above flotation. " *
    "Polar stereographic projection (70°N, 45°W)."

function plot_greenland(d, style; labels=nothing, maxtier=1, scale=1.85, ocean=:light, contours=true,
                         velcmap=CS_VEL)
    # Map window (km) fills the page width; beyond the BedMachine domain the
    # topography/bathymetry is ETOPO 2022 and glacier ice is from RGI 6.0.
    pad  = 70
    W    = A0_PORTRAIT[1] - 2pad
    H    = 2650.0
    ymap = (-3420.0, -600.0)
    kmpp = (ymap[2] - ymap[1])/H
    xc   = 120.0                                  # centre of Greenland
    xl   = (xc - W*kmpp/2, xc + W*kmpp/2)

    r   = d.r
    img = compose(style, r; srflim=(0, 3300), bedlim=(-1500, 3000), ocean, velcmap)

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
    if contours || style == :surface
        elevation_contours!(ax, xs, ys, r.zs[ix, iy], d.grounded[ix, iy], 500:500:3000)
    end
    maskoutline!(ax, xs, ys, r.ice[ix, iy]; σ=1.5, color=(:black, 0.75), linewidth=1.5)
    basin_divides!(ax, d.xb, d.yb, d.basin, d.ice_b; σ=1.2, color=(:black, 0.6), lw=2.5)
    regions = basin_labels!(ax, d.xb, d.yb, d.basin, d.basin_names, d.ice_b; kmpp, fontsize=44, color=(:gray20, 0.8))

    seaice = seaice_edges!(ax, "greenland", ("03", "09"); ocean)
    graticule!(ax, PROJ_GRL, 60:5:80, -80:10:0; latrange=(58, 84), color=(:gray20, 0.4), lw=1.0,
               inside=p -> xl[1] <= p[1] <= xl[2] && ymap[1] <= p[2] <= ymap[2])
    scalebar!(ax, 700.0, ymap[1] + 60, 400.0; fontsize=26, color=(ocean === :dark ? :white : :black))

    if labels !== nothing
        draw_labels!(ax, labels, PROJ_GRL; layout=:coastal, kmpp, maxtier, scale, obstacles=regions,
                     icemask=r.ice, groundedmask=d.grounded, x=r.x, y=r.y, offset=50.0, maxlead=300.0, tmax=120.0,
                     seacolor=seacolor(style, ocean), limits=(xl..., ymap...))
    end

    foot = fig[3, 1] = GridLayout()
    bars = foot[1, 1] = GridLayout()
    cb = (; vertical=false, width=500, height=28, labelsize=30, ticklabelsize=26)
    style_colorbar!(bars[1, 1], style; velcmap, cb...)
    style == :bed || ocean_colorbar!(bars[2, 1], ocean; cb...)
    rowgap!(bars, 16)
    contour_entry = contours || style == :surface ?
        [(LineElement(color=(:gray20, 0.7), linewidth=0.9scale), "500 m surface contours")] : []
    symbol_legend!(foot[1, 2]; scale, domes=false, labelsize=26, nbanks=2, rowgap=2, tellheight=true,
                   extra=vcat(contour_entry, seaice))
    numbers_box!(foot[1, 3], d.attrs; title="Greenland Ice Sheet in numbers", fontsize=26, valign=:top)
    qr_code!(foot[1, 4]; size=150, fontsize=22, valign=:top)
    colgap!(foot, 60)
    Label(fig[4, 1], CREDITS_GRL; fontsize=22, color=:gray30, word_wrap=true, width=W - 200, justification=:left)
    rowgap!(fig.layout, 30)
    return fig
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_poster("greenland", plot_greenland)
end
