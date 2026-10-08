# Shared tools for ice-sheet poster maps (projection, shading, labels).

using CairoMakie, NCDatasets, ColorSchemes, Colors, DelimitedFiles, Statistics, Unicode, Printf
import JSON
import Contour as CT

include("paths.jl")
include("projection.jl")

# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

function readvar(path, name)
    NCDataset(path) do ds
        Float64.(coalesce.(ds[name][:, :], NaN))
    end
end

readaxis(path, name) = NCDataset(ds -> Float64.(ds[name][:]), path)

"Separable Gaussian smoothing (sigma in grid cells), NaN-aware."
function gauss_smooth(z, σ)
    r = ceil(Int, 3σ)
    k = [exp(-0.5*(i/σ)^2) for i in -r:r]
    function pass(a, dim)
        b = similar(a); nx, ny = size(a)
        for j in 1:ny, i in 1:nx
            s = 0.0; w = 0.0
            for (q, kq) in zip(-r:r, k)
                ii, jj = dim == 1 ? (i+q, j) : (i, j+q)
                (1 <= ii <= nx && 1 <= jj <= ny) || continue
                v = a[ii, jj]; isnan(v) && continue
                s += kq*v; w += kq
            end
            b[i, j] = w > 0 ? s/w : NaN
        end
        return b
    end
    return pass(pass(z, 1), 2)
end

# ---------------------------------------------------------------------------
# Shading and colour compositing
# ---------------------------------------------------------------------------

"Lambertian hillshade anomaly relative to flat ground (0 = flat); dx in km, zfac = vertical exaggeration."
function hillshade(z, dx; az=315.0, alt=40.0, zfac=1.0)
    nx, ny = size(z)
    hs = fill(NaN, nx, ny)
    a = deg2rad(az); h = deg2rad(alt)
    lx, ly, lz = cos(h)*sin(a), cos(h)*cos(a), sin(h)
    d = dx*1e3
    for j in 1:ny, i in 1:nx
        i0, i1 = max(i-1, 1), min(i+1, nx)
        j0, j1 = max(j-1, 1), min(j+1, ny)
        gx = zfac*(z[i1, j] - z[i0, j]) / ((i1 - i0)*d)
        gy = zfac*(z[i, j1] - z[i, j0]) / ((j1 - j0)*d)
        (isnan(gx) || isnan(gy)) && continue
        n = sqrt(gx^2 + gy^2 + 1)
        hs[i, j] = clamp((-gx*lx - gy*ly + lz)/n, 0, 1) - lz
    end
    return hs
end

"Map a value through a colormap with limits."
function cmap(cs, v, lo, hi)
    isnan(v) && return RGBf(NaN, NaN, NaN)
    return RGBf(get(cs, clamp((v - lo)/(hi - lo), 0, 1)))
end

"Brighten/darken colour by hillshade anomaly (0 = neutral)."
function shade(c::RGBf, hs; strength=1.0)
    isnan(hs) && return c
    f = 1 + strength*hs
    return RGBf(clamp(c.r*f, 0, 1), clamp(c.g*f, 0, 1), clamp(c.b*f, 0, 1))
end

mix(a::RGBf, b::RGBf, t) = RGBf(a.r + t*(b.r - a.r), a.g + t*(b.g - a.g), a.b + t*(b.b - a.b))

const CS_OCEAN = cgrad([colorant"#0a1f3d", colorant"#163d6b", colorant"#2f6ea6", colorant"#7fb3d9"], [0, 0.45, 0.8, 1.0])
"Colour of sea names for a raster style and ocean variant."
seacolor(style, ocean) = style == :bed ? :gray20 : (ocean === :light ? colorant"#1d4e89" : :white)

# light variant: shallow shelf light blue, deep ocean fading to white
const CS_OCEAN_LIGHT = cgrad([colorant"#ffffff", colorant"#e8f1f8", colorant"#c6ddef", colorant"#97c4e6"], [0, 0.45, 0.85, 1.0])
const CS_ICE   = cgrad([colorant"#5f86ad", colorant"#8eadca", colorant"#bfd2e3", colorant"#e3ecf4", colorant"#fbfdff"], [0, 0.2, 0.45, 0.75, 1.0])
const CS_ROCK  = cgrad([colorant"#6e604c", colorant"#a08c6c", colorant"#d2c4a8"])
const CS_BED_LO = cgrad(reverse(ColorSchemes.oleron.colors[1:128]))   # 0 → deep
const CS_BED_HI = cgrad(ColorSchemes.oleron.colors[129:end])           # 0 → high

"Bed elevation colour with asymmetric ranges about sea level."
bedcolor(v, lo, hi) = v < 0 ? cmap(CS_BED_LO, -v, 0, -lo) : cmap(CS_BED_HI, v, 0, hi)

# ---------------------------------------------------------------------------
# Contours (Contour.jl; NaN-safe, unlike Makie's contour! with labels machinery)
# ---------------------------------------------------------------------------

"""
Contour line vertices for given levels, NaN-separated. `z` must be finite;
vertices whose nearest grid node is outside `mask` are dropped (line breaks).
"""
function contour_points(x, y, z, levels; mask=nothing)
    pts = Point2f[]
    dx = x[2] - x[1]; dy = y[2] - y[1]
    for lv in levels, ln in CT.lines(CT.contour(x, y, z, lv))
        xs, ys = CT.coordinates(ln)
        for (a, b) in zip(xs, ys)
            ok = true
            if mask !== nothing
                i = clamp(round(Int, (a - x[1])/dx) + 1, 1, length(x))
                j = clamp(round(Int, (b - y[1])/dy) + 1, 1, length(y))
                ok = mask[i, j]
            end
            push!(pts, ok ? Point2f(a, b) : Point2f(NaN, NaN))
        end
        push!(pts, Point2f(NaN, NaN))
    end
    return pts
end


"Surface-elevation contours over `mask`, from a surface smoothed by `σkm` (vertices, NaN-separated)."
elevation_contour_points(x, y, zs, mask, levels; σkm=2.0) =
    contour_points(x, y, gauss_smooth(zs, σkm/(x[2] - x[1])), levels; mask)

"Thin surface-elevation contours over `mask`, from a surface smoothed by `σkm`."
elevation_contours!(ax, x, y, zs, mask, levels; σkm=2.0, color=(:gray20, 0.7), linewidth=0.9) =
    lines!(ax, elevation_contour_points(x, y, zs, mask, levels; σkm); color, linewidth)

"Smoothed outline of a boolean mask (vertices, NaN-separated); `mask` limits where it is drawn."
maskoutline_points(x, y, m; σ=0.7, mask=nothing) = contour_points(x, y, gauss_smooth(Float64.(m), σ), [0.5]; mask)
maskoutline!(ax, x, y, m; σ=0.7, mask=nothing, kw...) = lines!(ax, maskoutline_points(x, y, m; σ, mask); kw...)

# ---------------------------------------------------------------------------
# Map furniture
# ---------------------------------------------------------------------------

"Draw graticule lines clipped to the axis limits."
function graticule!(ax, proj, lats, lons; latrange, lonres=0.25, color=(:white, 0.35), lw=0.6,
                    inside=p -> true)
    clip(xy) = [inside(p) ? p : Point2f(NaN, NaN) for p in xy]
    for φ in lats
        xy = [Point2f(project(proj, λ, φ)...) for λ in -180:lonres:180]
        lines!(ax, clip(xy); color=color, linewidth=lw, linestyle=:dash)
    end
    for λ in lons
        xy = [Point2f(project(proj, λ, φ)...) for φ in range(latrange..., length=400)]
        lines!(ax, clip(xy); color=color, linewidth=lw, linestyle=:dash)
    end
end

"Scale bar in km at lower-left corner position (x0,y0)."
function scalebar!(ax, x0, y0, len; segs=4, h=len/40, fontsize=14, color=:black)
    w = len/segs
    for k in 0:segs-1
        poly!(ax, Rect(x0 + k*w, y0, w, h); color=(isodd(k) ? :white : :black), strokecolor=:black, strokewidth=0.8)
    end
    text!(ax, x0, y0 + 1.8h; text="0", fontsize, color, align=(:center, :bottom))
    text!(ax, x0 + len, y0 + 1.8h; text="$(round(Int, len)) km", fontsize, color, align=(:center, :bottom))
end

"Interfaces between integer basin ids over `icemask` (vertices, NaN-separated)."
function basin_divide_points(x, y, basin, icemask; σ=1.0)
    pts = Point2f[]
    for id in sort(unique(filter(>(0), basin)))
        append!(pts, contour_points(x, y, gauss_smooth(Float64.(basin .== id), σ), [0.5]; mask=icemask))
    end
    return pts
end

basin_divides!(ax, x, y, basin, icemask; σ=1.0, color=(:black, 0.55), lw=1.0, linestyle=:solid) =
    lines!(ax, basin_divide_points(x, y, basin, icemask; σ); color, linewidth=lw, linestyle)

# ---------------------------------------------------------------------------
# Raster composition for the three styles
# ---------------------------------------------------------------------------

# Velocity colour maps (log10 speed, 0.3-3.5): slow ice fades into the hillshade
const VEL_CMAPS = Dict(
    # classic Rignot/Mouginot-style: beige-green-cyan-blue-purple-magenta
    :classic => cgrad([colorant"#f3e7c4", colorant"#9fd27a", colorant"#43b8b4",
                       colorant"#2a63b3", colorant"#7a2fbf", colorant"#e3238c"]),
    # cold to hot: slow ice in cool blues, fast outlets glowing orange/yellow
    :ember   => cgrad([colorant"#dcecf2", colorant"#8cc5d9", colorant"#3f84b8", colorant"#2b3a83",
                       colorant"#9b2f74", colorant"#e8622c", colorant"#ffd23f"]),
    # perceptually uniform, colour-blind safe (Crameri batlow)
    :batlow  => cgrad(ColorSchemes.batlow.colors),
    # perceptually uniform, single-hue-ish fire (Crameri lajolla, light to dark)
    :lajolla => cgrad(ColorSchemes.lajolla.colors),
)
const CS_VEL = VEL_CMAPS[:classic]

const C_SHELF = RGBf(colorant"#b9cbdb")

"""
Raster fields for `compose`: hillshades of surface and bed (dx in km) plus
masks derived from the prepared mask (0 ocean, 1 land, 2 grounded, 3 floating).
"""
function raster_fields(x, y, zs, zb, mask, u; dx, zfac_ice=40, zfac_bed=6)
    return (; x, y, dx, zs, zb, u,
              hs_s=hillshade(zs, dx; zfac=zfac_ice), hs_b=hillshade(zb, dx; zfac=zfac_bed),
              ice=mask .>= 2, ocean=mask .== 0, shelf=mask .== 3)
end

"""
Read a prepared poster grid (scripts/prepare.jl): the full-resolution grid in
\$CRYOMAPS_DATA (`full=true`), or its copy at every DRAFT_STRIDE-th node in the
repo. Returns the raster fields, the basin ids on the coarser basin grid with
matching coarse masks, and the global attributes (incl. the key numbers).
"""
function load_prepared(region; zfac_ice=40, zfac_bed=6, full=false)
    f = full ? prepared_file(region) : draft_file(region)
    isfile(f) || error("missing $f: run steps 0-1 (see README)")
    x = readaxis(f, "x"); y = readaxis(f, "y")
    xb = readaxis(f, "xb"); yb = readaxis(f, "yb")
    mask = NCDataset(ds -> Int.(ds["mask"][:, :]), f)
    r = raster_fields(x, y, readvar(f, "z_srf"), readvar(f, "z_bed"), mask, readvar(f, "u");
                      dx=x[2] - x[1], zfac_ice, zfac_bed)
    st = round(Int, (xb[2] - xb[1])/(x[2] - x[1]))       # basin grid nodes are every st-th node
    basin, basin_names = NCDataset(ds -> (Int.(ds["basin"][:, :]), split(ds["basin"].attrib["names"], ",")), f)
    return (; r, xb, yb, basin, basin_names, grounded_b=(mask[1:st:end, 1:st:end] .== 2), ice_b=(mask[1:st:end, 1:st:end] .>= 2),
              grounded=(mask .== 2), attrs=NCDataset(ds -> Dict(ds.attrib), f))
end

"""
RGB image for a given style from `raster_fields` output:
  :surface  – ice-surface elevation tint + hillshade, bathymetry, rock
  :bed      – bed topography everywhere (oleron), ice margin drawn separately
  :velocity – grey hillshaded surface with log-velocity overlay
`ocean` is :dark (deep navy) or :light (deep ocean fading to white).
"""
function compose(style, r; srflim=(0, 3300), bedlim=(-2000, 3000), oceanlim=(-4500, 0), ocean=:light,
                 velcmap=CS_VEL)
    cs_ocean, hs_ocean = ocean === :light ? (CS_OCEAN_LIGHT, 0.25) : (CS_OCEAN, 0.5)
    nx, ny = size(r.zs)
    img = Matrix{RGBf}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        ice, ocean = r.ice[i, j], r.ocean[i, j] && !r.ice[i, j]
        if style == :bed
            c = shade(bedcolor(r.zb[i, j], bedlim...), r.hs_b[i, j]; strength=0.9)
            ocean && (c = mix(c, RGBf(1, 1, 1), 0.25))
        elseif ocean
            c = shade(cmap(cs_ocean, r.zb[i, j], oceanlim...), r.hs_b[i, j]; strength=hs_ocean)
        elseif !ice
            c = shade(cmap(CS_ROCK, r.zs[i, j], 0, 2000), r.hs_b[i, j]; strength=1.0)
        elseif style == :surface && r.shelf[i, j]
            c = shade(C_SHELF, r.hs_s[i, j]; strength=0.6)
        elseif style == :surface
            c = shade(cmap(CS_ICE, r.zs[i, j], srflim...), r.hs_s[i, j]; strength=0.9)
        else  # :velocity
            c = shade(RGBf(0.86, 0.86, 0.86), r.hs_s[i, j]; strength=0.8)
            u = r.u[i, j]
            if !isnan(u) && u > 0
                lu = log10(max(u, 1.0))
                α = clamp((lu - 0.3)/0.9, 0, 1)
                c = mix(c, cmap(velcmap, lu, 0.3, 3.5), 0.85α)
            end
        end
        img[i, j] = c
    end
    return img
end

"Label point of each basin over ice: the in-basin point nearest its centroid. Returns (name, (x, y)) pairs."
function basin_label_points(x, y, basin, names, icemask)
    out = Tuple{String, Tuple{Float64, Float64}}[]
    for id in sort(unique(filter(>(0), basin[icemask])))
        I = findall((basin .== id) .& icemask)
        length(I) < 30 && continue
        cx = mean(x[c[1]] for c in I); cy = mean(y[c[2]] for c in I)
        k = argmin([hypot(x[c[1]] - cx, y[c[2]] - cy) for c in I])
        push!(out, (String(names[id]), (x[I[k][1]], y[I[k][2]])))
    end
    return out
end

"""
Label each basin over ice with its name (see basin_label_points). Returns the
label boxes (km) so place-name labels can avoid them.
"""
function basin_labels!(ax, x, y, basin, names, icemask; kmpp, fontsize=16, color=:gray20)
    boxes = Tuple[]
    for (name, p) in basin_label_points(x, y, basin, names, icemask)
        halotext!(ax, p[1], p[2]; text=name, fontsize, color, font=:bold, align=(:center, :center))
        push!(boxes, boxat(p, textbox(name, fontsize, kmpp; font=:bold), :center))
    end
    return boxes
end

# ---------------------------------------------------------------------------
# Poster furniture (sizes in pt; figures are saved with 1 unit = 1 pt)
# ---------------------------------------------------------------------------

const A0_PORTRAIT  = (2384, 3370)
const A0_LANDSCAPE = (3370, 2384)

"Colorbar matching the raster style."
function style_colorbar!(pos, style; bedlim=(-1.5, 3.0), srflim=(0, 3.3), velcmap=CS_VEL, kw...)
    if style == :velocity
        Colorbar(pos; colormap=velcmap, limits=(0.3, 3.5), ticks=(0:3, ["1", "10", "100", "1000"]),
                 label="Surface ice velocity [m/yr]", kw...)
    elseif style == :surface
        Colorbar(pos; colormap=CS_ICE, limits=srflim, label="Ice-surface elevation [km]", kw...)
    else
        Colorbar(pos; colormap=cgrad(vcat(reverse(CS_BED_LO.colors.colors[1:4:end]), CS_BED_HI.colors.colors[1:2:end])),
                 limits=bedlim, label="Bed elevation [km]", kw...)
    end
end

"Ocean colour scale (depth in m) for the light or dark ocean."
function ocean_colorbar!(pos, ocean; oceanlim=(-4500, 0), kw...)
    cs = ocean === :light ? CS_OCEAN_LIGHT : CS_OCEAN
    Colorbar(pos; colormap=reverse(cs), limits=(0, -oceanlim[1]), ticks=0:1000:-oceanlim[1],
             label="Ocean depth [m]", kw...)
end

# ---------------------------------------------------------------------------
# Sea-ice edge, key numbers, QR code
# ---------------------------------------------------------------------------

const MONTHS = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
                "October", "November", "December"]

"Line vertices (km, NaN-separated) of a GeoJSON file of (multi)line strings in metres."
function geojson_lines(path)
    pts = Point2f[]
    for f in JSON.parsefile(path)["features"]
        g = f["geometry"]; g === nothing && continue
        parts = g["type"] == "LineString" ? [g["coordinates"]] : g["coordinates"]
        for ln in parts
            append!(pts, [Point2f(c[1]/1e3, c[2]/1e3) for c in ln]); push!(pts, Point2f(NaN, NaN))
        end
    end
    return pts
end

"""
Median sea-ice edges (1981-2010) for `months` ("MM", winter maximum first),
dashed for the maximum and dotted for the minimum. Returns legend entries.
"""
function seaice_edges!(ax, region, months; ocean=:light, linewidth=2.5)
    col = ocean === :light ? (colorant"#1d4e89", 0.8) : (colorant"#7cc8f5", 0.9)    # visible on map and legend
    entries = []
    for (mm, ls) in zip(months, (:dash, :dot))
        f = joinpath(REPO_PREP_DIR, "seaice_$(region)_$(mm).geojson")
        lines!(ax, geojson_lines(f); color=col, linewidth, linestyle=ls)
        push!(entries, (LineElement(; color=col, linewidth, linestyle=ls),
                        "Sea-ice edge, $(MONTHS[parse(Int, mm)])"))
    end
    return entries
end

"Key numbers of the ice sheet (from prepare.jl) as a two-column table."
function numbers_box!(pos, attrs; title="Key numbers", fontsize=28, kw...)
    a(k) = Float64(attrs[k])
    rows = ["Ice area"             => @sprintf("%.2f million km²", a("ice_area_km2")/1e6),
            "   of which floating" => @sprintf("%.2f million km²", a("floating_area_km2")/1e6),
            "Ice volume"           => @sprintf("%.2f million km³", a("ice_volume_km3")/1e6),
            "Maximum thickness"    => @sprintf("%.0f m", a("max_thickness_m")),
            "Sea-level equivalent" => @sprintf("%.1f m", a("sea_level_equivalent_m"))]
    a("floating_area_km2") < 1e4 && deleteat!(rows, 2)
    g = GridLayout(pos; kw...)
    Label(g[1, 1:2], title; fontsize=1.15fontsize, font=:bold, halign=:left)
    for (i, (k, v)) in enumerate(rows)
        Label(g[i+1, 1], k; fontsize, halign=:left, color=:gray25)
        Label(g[i+1, 2], v; fontsize, halign=:right, font=:bold)
    end
    colgap!(g, 30); rowgap!(g, 4)
    return g
end

"QR code of the website (data/qr_site.txt, see scripts/make_qr.jl) with a caption."
function qr_code!(pos; size=170, fontsize=24, kw...)
    ln = readlines(joinpath(ROOT, "data", "qr_site.txt"))
    url = strip(ln[1][2:end])
    m = permutedims(reduce(hcat, [[c == '1' for c in l] for l in ln[2:end]]))      # rows top to bottom
    n = Base.size(m, 1)
    g = GridLayout(pos; kw...)
    ax = Axis(g[1, 1]; width=size, height=size, aspect=1, limits=(-2, n + 2, -2, n + 2))
    hidedecorations!(ax); hidespines!(ax)
    heatmap!(ax, 0.5:1:n, 0.5:1:n, reverse(permutedims(m), dims=2); colormap=[:white, :black], colorrange=(0, 1))
    Label(g[2, 1], "Interactive map:\n" * replace(url, r"^https://" => "", r"/$" => ""); fontsize, color=:gray25)
    rowgap!(g, 6)
    return g
end

"Legend of map symbols. `extra` adds (element, label) pairs."
function symbol_legend!(pos; scale=1.5, domes=true, extra=[], rowgap=6, kw...)
    ms = 8*scale*1.6
    els = Any[MarkerElement(marker=:diamond, color=:darkred, strokecolor=:white, markersize=ms),
              MarkerElement(marker=:circle, color=:black, strokecolor=:white, markersize=ms)]
    labs = Any["Ice core", "Station / settlement"]
    if domes
        push!(els, MarkerElement(marker=:utriangle, color=:black, strokecolor=:white, markersize=ms)); push!(labs, "Dome")
    end
    append!(els, [LineElement(color=(:black, 0.6), linewidth=2.5), LineElement(color=(:black, 0.75), linewidth=1.5)])
    append!(labs, ["Drainage divide", "Ice margin"])
    for (e, l) in extra
        push!(els, e); push!(labs, l)
    end
    Legend(pos, els, labs; framevisible=false, rowgap, patchsize=(40, 20), kw...)
end

"""
Save the poster as a small PNG for sharing (`small_px` pixels on the long
edge) and, with `full`, also as a PDF at true size (1 unit = 1 pt) and a PNG
preview at `dpi_png`.
"""
function save_poster(fig, out; full=false, dpi_png=100, small_px=1600)
    mkpath(dirname(out))
    if full
        save(out*".pdf", fig; pt_per_unit=1)
        save(out*".png", fig; px_per_unit=dpi_png/72)
    end
    save(out*"_small.png", fig; px_per_unit=small_px/maximum(size(fig.scene)))
    println("saved ", out, full ? ".{pdf,png} and _small.png" : "_small.png")
end

"""
Command-line options:
- a raster style: velocity (default) | surface | bed
- `dark`: dark ocean; `nocontours`: no surface contours
- `cmap=<name>`: velocity colour map (see VEL_CMAPS)
- `tier=2`: also show the tier-2 names of the label CSV
- `add=<name>;<name>...`: add names from the label CSV or the gazetteer
- `list` or `list=<text>`: print the available names (matching <text>) and exit
- `full`: full-resolution grid (\$CRYOMAPS_DATA) and PDF + PNG output; without
  it, only the small PNG is made from the coarse grid in data/prepared.
"""
function parse_args(args)
    val(key) = (i = findfirst(startswith(key*"="), args); i === nothing ? nothing : split(args[i], "="; limit=2)[2])
    known = ("velocity", "surface", "bed", "dark", "nocontours", "full", "list")
    for a in args
        a in known || any(startswith(a, k*"=") for k in ("cmap", "tier", "add", "list")) || error("unknown option $a")
    end
    i = findfirst(in(("velocity", "surface", "bed")), args)
    style = Symbol(i === nothing ? "velocity" : args[i])
    ocean = "dark" in args ? :dark : :light
    contours = !("nocontours" in args)
    cmapname = Symbol(something(val("cmap"), "classic"))
    haskey(VEL_CMAPS, cmapname) || error("unknown cmap $cmapname; options: $(keys(VEL_CMAPS))")
    tier = parse(Int, something(val("tier"), "1"))
    tier in (1, 2) || error("tier must be 1 or 2; add gazetteer names with add=<name>")
    add = filter(!isempty, strip.(split(something(val("add"), ""), ";")))
    list = "list" in args ? "" : val("list")
    tag = join(filter(!isempty, [ocean === :dark ? "dark" : "", contours ? "" : "nocontours",
                                 cmapname === :classic ? "" : String(cmapname), tier == 1 ? "" : "tier$tier",
                                 isempty(add) ? "" : "custom"]), "_")
    return (; style, ocean, contours, velcmap=VEL_CMAPS[cmapname], tier, add, list, full="full" in args,
              tag=isempty(tag) ? "" : "_"*tag)
end

"Output path (without extension): the default poster in plots/, all others in plots/variants/."
function poster_path(region, o)
    name = "$(region)_A0_$(o.style)$(o.tag)"
    return o.style === :velocity && isempty(o.tag) ? joinpath(ROOT, "plots", name) : joinpath(ROOT, "plots", "variants", name)
end

"""
Command-line entry point of the poster scripts: `plotfun(d, style; labels,
maxtier, ocean, contours, velcmap)` draws the poster of `region`.
"""
function main_poster(region, plotfun, args=ARGS)
    o = parse_args(args)
    labs = read_labels(joinpath(ROOT, "data", "labels_$(region).csv"))
    gaz  = read_labels(joinpath(REPO_PREP_DIR, "gazetteer_$(region).csv"))
    if o.list !== nothing
        list_names(vcat(labs, gaz), o.list); return
    end
    labs = add_names(labs, gaz, o.add)
    d = load_prepared(region; o.full)
    fig = plotfun(d, o.style; labels=labs, maxtier=o.tier, o.ocean, o.contours, o.velcmap)
    save_poster(fig, poster_path(region, o); o.full)
end

include("labels.jl")
