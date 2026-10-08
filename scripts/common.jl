# Shared tools for ice-sheet poster maps (projection, shading, labels).

using CairoMakie, NCDatasets, ColorSchemes, Colors, DelimitedFiles, Statistics
import Contour as CT

include("paths.jl")

# ---------------------------------------------------------------------------
# Projection: ellipsoidal polar stereographic (Snyder 1987), output in km
# ---------------------------------------------------------------------------

struct PolarStereo
    lon0::Float64     # straight vertical longitude from pole [deg]
    lat_ts::Float64   # latitude of true scale [deg] (sign gives hemisphere)
end

const WGS84_A = 6378137.0
const WGS84_E = 0.0818191908426

function _tfun(φ, e)
    return tan(π/4 - φ/2) / ((1 - e*sin(φ)) / (1 + e*sin(φ)))^(e/2)
end

function project(p::PolarStereo, lon, lat)
    s  = sign(p.lat_ts)                  # +1 north, -1 south
    φ  = deg2rad(s*lat)
    φc = deg2rad(s*p.lat_ts)
    λ  = deg2rad(s*(lon - p.lon0))
    e  = WGS84_E
    mc = cos(φc) / sqrt(1 - e^2*sin(φc)^2)
    ρ  = WGS84_A * mc * _tfun(φ, e) / _tfun(φc, e)
    x  =  ρ*sin(λ)
    y  = -ρ*cos(λ)
    return (s*x/1e3, s*y/1e3)
end

const PROJ_GRL = PolarStereo(-45.0, 70.0)
const PROJ_ANT = PolarStereo(0.0, -71.0)

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

contourlines!(ax, x, y, z, levels; mask=nothing, kw...) = lines!(ax, contour_points(x, y, z, levels; mask); kw...)

"Smoothed outline of a boolean mask."
maskoutline!(ax, x, y, m; σ=0.7, kw...) = contourlines!(ax, x, y, gauss_smooth(Float64.(m), σ), [0.5]; kw...)

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

"Contour lines of the interfaces between integer basin ids, drawn only over ice."
function basin_divides!(ax, x, y, basin, icemask; σ=1.0, color=(:black, 0.55), lw=1.0, linestyle=:solid)
    ids = sort(unique(filter(>(0), basin)))
    for id in ids
        ind = Float64.(basin .== id)
        ind = gauss_smooth(ind, σ)
        contourlines!(ax, x, y, ind, [0.5]; mask=icemask, color=color, linewidth=lw, linestyle=linestyle)
    end
end

# ---------------------------------------------------------------------------
# Raster composition for the three styles
# ---------------------------------------------------------------------------

const CS_VEL = cgrad([colorant"#f3e7c4", colorant"#9fd27a", colorant"#43b8b4",
                      colorant"#2a63b3", colorant"#7a2fbf", colorant"#e3238c"])

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
Read a prepared poster grid (scripts/prepare.jl). Returns the raster fields and,
on the basin grid, the basin ids with a matching coarse ice mask.
"""
function load_prepared(region; zfac_ice, zfac_bed=6)
    f = prepared_file(region)
    isfile(f) || error("missing $f: run steps 0-1 (see README)")
    x = readaxis(f, "x"); y = readaxis(f, "y"); xb = readaxis(f, "xb"); yb = readaxis(f, "yb")
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
"""
function compose(style, r; srflim=(0, 3300), bedlim=(-2000, 3000), oceanlim=(-4500, 0))
    nx, ny = size(r.zs)
    img = Matrix{RGBf}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        ice, ocean = r.ice[i, j], r.ocean[i, j] && !r.ice[i, j]
        if style == :bed
            c = shade(bedcolor(r.zb[i, j], bedlim...), r.hs_b[i, j]; strength=0.9)
            ocean && (c = mix(c, RGBf(1, 1, 1), 0.25))
        elseif ocean
            c = shade(cmap(CS_OCEAN, r.zb[i, j], oceanlim...), r.hs_b[i, j]; strength=0.5)
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
                c = mix(c, cmap(CS_VEL, lu, 0.3, 3.5), 0.85α)
            end
        end
        img[i, j] = c
    end
    return img
end

"""
Label each basin over ice with its name, at the in-basin point nearest its
centroid. Returns the label boxes (km) so place-name labels can avoid them.
"""
function basin_labels!(ax, x, y, basin, names, icemask; kmpp, fontsize=16, color=:gray20)
    boxes = Tuple[]
    for id in sort(unique(filter(>(0), basin[icemask])))
        I = findall((basin .== id) .& icemask)
        length(I) < 30 && continue
        cx = mean(x[c[1]] for c in I); cy = mean(y[c[2]] for c in I)
        k = argmin([hypot(x[c[1]] - cx, y[c[2]] - cy) for c in I])
        p = (x[I[k][1]], y[I[k][2]])
        halotext!(ax, p[1], p[2]; text=names[id], fontsize, color, font=:bold, align=(:center, :center))
        push!(boxes, boxat(p, textbox(names[id], fontsize, kmpp; font=:bold), :center))
    end
    return boxes
end

# ---------------------------------------------------------------------------
# Poster furniture (sizes in pt; figures are saved with 1 unit = 1 pt)
# ---------------------------------------------------------------------------

const A0_PORTRAIT  = (2384, 3370)
const A0_LANDSCAPE = (3370, 2384)

"Colorbar matching the raster style."
function style_colorbar!(pos, style; bedlim=(-1.5, 3.0), srflim=(0, 3.3), kw...)
    if style == :velocity
        Colorbar(pos; colormap=CS_VEL, limits=(0.3, 3.5), ticks=(0:3, ["1", "10", "100", "1000"]),
                 label="Surface ice velocity [m/yr]", kw...)
    elseif style == :surface
        Colorbar(pos; colormap=CS_ICE, limits=srflim, label="Ice-surface elevation [km]", kw...)
    else
        Colorbar(pos; colormap=cgrad(vcat(reverse(CS_BED_LO.colors.colors[1:4:end]), CS_BED_HI.colors.colors[1:2:end])),
                 limits=bedlim, label="Bed elevation [km]", kw...)
    end
end

"Legend of map symbols. `extra` adds (element, label) pairs."
function symbol_legend!(pos; scale=1.5, domes=true, extra=[], kw...)
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
    Legend(pos, els, labs; framevisible=false, rowgap=6, patchsize=(40, 20), kw...)
end

"Save the poster PDF at true size (1 unit = 1 pt) and a PNG preview."
function save_poster(fig, out; dpi_png=100)
    save(out*".pdf", fig; pt_per_unit=1)
    save(out*".png", fig; px_per_unit=dpi_png/72)
    println("saved ", out, ".{pdf,png}")
end

include("labels.jl")
