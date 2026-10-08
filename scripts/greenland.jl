# Greenland A0 portrait poster map.
# Usage: julia --project=. scripts/greenland.jl [velocity|surface|bed ...]
# Requires data/GRL-8KM-EXT_ETOPO2022.nc (see README).

include("common.jl")

const GDIR = joinpath(ICE_DATA, "Greenland", "GRL-8KM")
const DX   = 8.0

const CREDITS_GRL = "Data: bed and surface topography from BedMachine Greenland v3 (Morlighem et al., 2017); " *
    "surface ice velocity from Joughin et al. (2018); drainage basins (numbered) from Zwally et al. (2012); " *
    "outside the BedMachine domain, topography from ETOPO 2022 (NOAA NCEI, 2022) and glacier outlines " *
    "from RGI 6.0 (RGI Consortium, 2017); " *
    "place names from GeoNames. All fields on an 8 km polar stereographic grid (70°N, 45°W)."

const ETOPO_EXT = joinpath(ROOT, "data", "GRL-8KM-EXT_ETOPO2022.nc")
const RGI_EXT   = joinpath(ROOT, "data", "GRL-8KM-EXT_RGI60.nc")

"""
Load GRL-8KM fields and embed them in the extended (padded) 8 km grid made by
prepare_etopo.jl. Outside the BedMachine domain, ETOPO 2022 gives bed/surface
elevation and RGI 6.0 (ice_frac > 0.5) marks glacier ice; there is no velocity
or basin information there.
"""
function load_greenland(; f=8)
    for f in (ETOPO_EXT, RGI_EXT)
        isfile(f) || error("missing $(f): run the data preparation steps in README.md")
    end
    fn = joinpath(GDIR, "GRL-8KM_TOPO-M17.nc")
    xg = readaxis(fn, "xc"); yg = readaxis(fn, "yc")
    x  = readaxis(ETOPO_EXT, "xc"); y = readaxis(ETOPO_EXT, "yc")
    zx = readvar(ETOPO_EXT, "z")
    zx[isnan.(zx)] .= minimum(filter(!isnan, zx))      # empty cells lie outside the poster window

    # extended fields from ETOPO, then overwrite the BedMachine domain
    zb = copy(zx); zs = max.(zx, 0.0); ocean = zx .< 0
    H  = Float64.(readvar(RGI_EXT, "ice_frac") .> 0.5)   # nominal thickness: only the ice mask is used
    u = fill(NaN, size(zx)); basin = zeros(size(zx))
    I = round(Int, (xg[1] - x[1])/DX) .+ (1:length(xg))
    J = round(Int, (yg[1] - y[1])/DX) .+ (1:length(yg))
    zb[I, J] = readvar(fn, "z_bed"); zs[I, J] = readvar(fn, "z_srf"); H[I, J] = readvar(fn, "H_ice")
    ocean[I, J] = readvar(fn, "mask") .== 0
    u[I, J] = readvar(joinpath(GDIR, "GRL-8KM_VEL-J18.nc"), "uxy_srf")
    basin[I, J] = readvar(joinpath(GDIR, "GRL-8KM_BASINS-nasa.nc"), "basin")

    ice = H .> 0
    r = raster_fields(x, y, zs, zb, H, ocean, falses(size(H)), u; dx=DX, f, zfac_ice=60, zfac_bed=6)
    return (; x, y, ice, basin, r)
end

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
    Label(head[2, 1], "Surface ice velocity, drainage basins and place names"; fontsize=44, color=:gray30)
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
    basin_divides!(ax, d.x, d.y, d.basin, d.ice; σ=1.2, color=(:black, 0.6), lw=2.5)
    basin_ids!(ax, d.x, d.y, d.basin, d.ice; fontsize=34)

    graticule!(ax, PROJ_GRL, 60:5:80, -80:10:0; latrange=(58, 84), color=(:gray20, 0.4), lw=1.0,
               inside=p -> xl[1] <= p[1] <= xl[2] && ymap[1] <= p[2] <= ymap[2])
    scalebar!(ax, 700.0, ymap[1] + 60, 400.0; fontsize=26, color=:white)

    if labels !== nothing
        draw_labels!(ax, labels, PROJ_GRL; layout=:coastal, kmpp, maxtier, scale,
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
