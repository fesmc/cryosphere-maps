# Paleo A0 landscape poster maps: the ice sheets at a past time slice, from PaleoMIST 1.0
# changes applied to present-day topography (grids from scripts/prepare_paleo.jl).
# Usage: julia --project=. scripts/paleo.jl nh|antarctica [surface|bed] [dark] [nocontours] [time=20] [pdf]
# Options as for the present-day posters (see parse_args in common.jl and the README);
# `time` is the time slice in ka (20 by default).

include("common.jl")

const PALEO_ROWS = [("Ice area [10⁶ km²]", "ice_area_km2", 1e6, "%.2f"),
                    ("Ice volume [10⁶ km³]", "ice_volume_km3", 1e6, "%.2f"),
                    ("Max. thickness [m]", "max_thickness_m", 1, "%.0f"),
                    ("Sea-level equiv. [m]", "sea_level_equivalent_m", 1, "%.1f"),
                    ("   more than today [m]", "sle_change_m", 1, "%.1f")]

const C_COAST = (colorant"#1d4e89", 0.8)       # present-day coastline
const C_PDICE = (:gray15, 0.8)                 # present-day ice margin

# Names of time slices shown in the subtitle
const PERIODS = Dict(20.0 => "Last Glacial Maximum")

"Present-day base data of the credits."
const PALEO_BASE = Dict(
    "nh" => "present-day topography and bathymetry from ETOPO 2022 (NOAA NCEI, 2022)",
    "antarctica" => "present-day bed and surface topography from BedMachine Antarctica v4 (Morlighem et al., 2020)")

const PALEO = Dict(
    "nh" => (proj=PROJ_GRL, xmap=(-5300.0, 4500.0), ymap=(-5700.0, 2500.0), title="Northern\nHemisphere\nIce Sheets",
             lims=(srflim=(0, 4100), bedlim=(-2500, 3000), oceanlim=(-4500, 0)),
             lats=30:10:80, lons=-180:30:150, latrange=(25, 88), labels="labels_paleo_nh.csv",
             cols=[("All", ""), ("N.Am.", "namerica_"), ("Grl.", "greenland_"), ("Eur.", "eurasia_")],
             note="N.Am.(erica), Gr(eenl)and, Eur.(asia): approximate regions;\nAll includes Iceland",
             projection="Polar stereographic projection (70°N, 45°W).", bottom=false),
    "antarctica" => (proj=PROJ_ANT, xmap=(-3040.0, 3040.0), ymap=(-2600.0, 2500.0), title="Antarctic\nIce Sheet",
             lims=STYLE_LIMITS["antarctica"],
             lats=-80:10:-60, lons=-180:30:150, latrange=(-88, -55), labels="labels_antarctica.csv",
             cols=[("All", ""), ("East", "east_"), ("West", "west_"), ("Pen.", "peninsula_")],
             note="East, West, Pen.(insula): IMBIE 2 regions, with ice\ngiven to the nearest region",
             projection="Polar stereographic projection (71°S, 0°E).", bottom=true),
)

"Years before present, with a thin space between thousands: 20 -> \"20 000\"."
years_bp(t) = (n = string(round(Int, 1000t)); length(n) > 3 ? n[1:end-3] * " " * n[end-2:end] : n)

credits(region, t) = "Data: ice sheets $(years_bp(t)) years ago from PaleoMIST 1.0 (Gowan et al., 2021): " *
    "the changes in ice thickness and in bed elevation relative to sea level (Earth deformation and sea-level " *
    "change) since then, added to $(PALEO_BASE[region]). PaleoMIST reconstructs grounded ice only. Key numbers " *
    "computed on the map grid; sea-level equivalent of the ice above flotation. " * PALEO[region].projection

"Read a paleo grid (scripts/prepare_paleo.jl) as raster fields for `compose`, with the present-day mask."
function load_paleo(region, t; zfac_ice=40, zfac_bed=6)
    f = paleo_file(region, t)
    isfile(f) || error("missing $f: run scripts/prepare_paleo.jl time=$t $region (see README)")
    x = readaxis(f, "x"); y = readaxis(f, "y")
    mask, pd_mask = NCDataset(ds -> (Int.(ds["mask"][:, :]), Int.(ds["pd_mask"][:, :])), f)
    r = raster_fields(x, y, readvar(f, "z_srf"), readvar(f, "z_bed"), mask, fill(NaN, size(mask));
                      dx=x[2] - x[1], zfac_ice, zfac_bed)
    return (; r, grounded=(mask .== 2), pd_mask, attrs=NCDataset(ds -> Dict(ds.attrib), f))
end

function plot_paleo(d, region, t, style; labels=nothing, scale=1.85, ocean=:light, contours=true)
    P = PALEO[region]
    xmap, ymap = P.xmap, P.ymap
    pad  = 50; panel = 520.0
    H    = A0_LANDSCAPE[2] - 2pad
    kmpp = (ymap[2] - ymap[1])/H
    W    = (xmap[2] - xmap[1])/kmpp

    r   = d.r
    img = compose(style, r; P.lims..., ocean)

    fig = Figure(size=A0_LANDSCAPE, figure_padding=pad, backgroundcolor=:white, fontsize=24)
    ax  = Axis(fig[1, 1]; width=W, height=H, limits=(xmap, ymap), backgroundcolor=:white)
    hidedecorations!(ax); hidespines!(ax)

    ix = findall(xi -> xmap[1] <= xi <= xmap[2], r.x); iy = findall(yi -> ymap[1] <= yi <= ymap[2], r.y)
    h = r.dx/2
    image!(ax, (r.x[ix[1]] - h, r.x[ix[end]] + h), (r.y[iy[1]] - h, r.y[iy[end]] + h), img[ix, iy]; interpolate=true)
    lines!(ax, Rect(xmap[1], ymap[1], xmap[2] - xmap[1], ymap[2] - ymap[1]); color=:black, linewidth=1.5)
    xs = r.x[ix]; ys = r.y[iy]

    if contours || style == :surface
        elevation_contours!(ax, xs, ys, r.zs[ix, iy], d.grounded[ix, iy], 500:500:4000; σkm=8.0)
    end
    pd = d.pd_mask[ix, iy]
    maskoutline!(ax, xs, ys, pd .>= 1; σ=0.7, color=C_COAST, linewidth=1.2, linestyle=:dash)
    maskoutline!(ax, xs, ys, pd .>= 2; σ=0.7, color=C_PDICE, linewidth=1.5, linestyle=:dot)
    maskoutline!(ax, xs, ys, r.ice[ix, iy]; σ=1.0, color=(:black, 0.75), linewidth=1.5)

    graticule!(ax, P.proj, P.lats, P.lons; latrange=P.latrange, color=(:gray20, 0.4), lw=1.0,
               inside=p -> xmap[1] <= p[1] <= xmap[2] && ymap[1] <= p[2] <= ymap[2])
    # scale bar and website in the lower or (NH: over the Pacific) upper left corner
    scalebar!(ax, xmap[1] + 70kmpp, P.bottom ? ymap[1] + 90kmpp : ymap[2] - 130kmpp, 1000.0; fontsize=26,
              color=(ocean === :dark ? :white : :black))
    url = site_url!(ax, (xmap..., ymap...), kmpp; corner=:left, P.bottom, ocean)

    if labels !== nothing
        draw_labels!(ax, labels, P.proj; layout=:coastal, kmpp, maxtier=1, scale, types=("region", "sea"), obstacles=[url],
                     icemask=r.ice, groundedmask=d.grounded, x=r.x, y=r.y, limits=(xmap..., ymap...),
                     seacolor=seacolor(style, ocean))
    end

    side = fig[1, 2] = GridLayout(width=panel)
    top  = side[1, 1] = GridLayout(valign=:top, tellheight=false)
    Label(top[1, 1], P.title; fontsize=96, font=:bold, halign=:left, justification=:left, lineheight=0.95)
    period = get(PERIODS, Float64(t), "")
    what = style == :bed ? "Bed topography" : "Ice-surface elevation"
    Label(top[2, 1], (isempty(period) ? "" : period * "\n") * "$(years_bp(t)) years ago\n$what";
          fontsize=36, color=:gray30, halign=:left, justification=:left)
    cb = (; vertical=false, width=panel - 40, height=28, labelsize=30, ticklabelsize=26, halign=:left)
    style_colorbar!(top[3, 1], style, r; P.lims.srflim, P.lims.bedlim, cb...)
    style == :bed || ocean_colorbar!(top[4, 1], ocean, r; P.lims.oceanlim, cb...)
    leg = [(LineElement(color=(:black, 0.75), linewidth=1.5), "Ice margin $(years_bp(t)) years ago"),
           (LineElement(color=C_PDICE, linewidth=1.5, linestyle=:dot), "Present-day ice margin"),
           (LineElement(color=C_COAST, linewidth=1.2, linestyle=:dash), "Present-day coastline")]
    (contours || style == :surface) &&
        push!(leg, (LineElement(color=(:gray20, 0.7), linewidth=0.9scale), "500 m surface contours"))
    Legend(top[5, 1], first.(leg), last.(leg); framevisible=false, rowgap=2, patchsize=(40, 20), labelsize=26,
           halign=:left, tellheight=true)
    numbers_table!(top[6, 1], d.attrs, P.cols; rows=PALEO_ROWS, title="$(years_bp(t)) years ago in numbers",
                   fontsize=24, halign=:left, note=P.note)
    rowgap!(top, 1, 20); rowgap!(top, 2, 70); rowgap!(top, 3, 20); rowgap!(top, 4, 60); rowgap!(top, 5, 60)
    bottom = side[2, 1] = GridLayout(valign=:bottom, tellheight=false)
    qr_code!(bottom[1, 1]; size=150, fontsize=22, halign=:left)
    Label(bottom[2, 1], credits(region, t); fontsize=22, color=:gray30, word_wrap=true, width=panel,
          justification=:left, halign=:left)
    rowgap!(bottom, 30)
    colgap!(fig.layout, 60)
    return fig
end

"""
Command-line entry point: the region (nh or antarctica), `time=<ka>` and the
options of the present-day posters that apply (surface or bed, dark, nocontours, pdf).
"""
function main_paleo(args=ARGS)
    reg = filter(in(keys(PALEO)), args)
    length(reg) == 1 || error("give one region: $(join(keys(PALEO), " or "))")
    region = only(reg)
    i = findfirst(startswith("time="), args)
    t = i === nothing ? 20.0 : parse(Float64, split(args[i], "=")[2])
    rest = filter(a -> a != region && !startswith(a, "time="), args)
    any(in(("surface", "bed")), rest) || push!(rest, "surface")
    o = parse_args(rest)
    o.style === :velocity && error("no velocity at past time slices; styles: surface, bed")
    o.full && error("`full`: the paleo grids exist at one resolution only (use `pdf` for PDF output)")
    labs = read_labels(joinpath(ROOT, "data", PALEO[region].labels))
    fig = plot_paleo(load_paleo(region, t), region, t, o.style; labels=labs, o.ocean, o.contours)
    save_poster(fig, joinpath(ROOT, "plots", "paleo", "$(region)_$(katag(t))_A0_$(o.style)$(o.tag)"); o.pdf)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_paleo()
end
