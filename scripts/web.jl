# Step 3: assets of the interactive web maps (Quarto site in site/).
#
# Writes to site/assets/:
#   <region>/tiles/<style>/<z>/<x>/<y>.jpg   raster tile pyramids of the map styles
#   <region>/*.geojson                       contours, ice margin, divides, names, sea-ice edges
#   <region>/config.json                     tile grid, projection, colour bars, key numbers
#   img/, posters/                           the small poster PNGs; PDFs and full PNGs
# and site/_gallery*.md, the poster galleries of the home, paleo and historical pages.
#
# Usage: julia --project=. scripts/web.jl [full] [greenland] [antarctica] [gallery]
# `full` uses the full-resolution grids in $CRYOMAPS_DATA (on albedo:
# jobs/web.sh); without it the coarse grids in data/prepared are used, which is
# enough to preview the site locally. `gallery` only updates the poster gallery.

include("common.jl")
using JpegTurbo
import Colors: N0f8

const SITE   = joinpath(ROOT, "site")
const ASSETS = joinpath(SITE, "assets")
const TILE   = 256

const PROJ4 = Dict(
    "greenland"  => "+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs",
    "antarctica" => "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
const PROJS = Dict("greenland" => PROJ_GRL, "antarctica" => PROJ_ANT)

# map styles of the web maps: (id, name, compose style, velocity colour map)
const WEB_STYLES = [("velocity", "Ice velocity", :velocity, :classic),
                    ("velocity_ember", "Ice velocity (ember)", :velocity, :ember),
                    ("surface", "Surface elevation", :surface, :classic),
                    ("bed", "Bed elevation", :bed, :classic)]
const SEAICE_MONTHS = Dict("greenland" => ("03", "09"), "antarctica" => ("09", "02"))

# ---------------------------------------------------------------------------
# Raster tiles
# ---------------------------------------------------------------------------

"Image (x, y) with ascending y as a row-major RGB{N0f8} matrix (top row first)."
rows_top_down(img) = [RGB{N0f8}(clamp(c.r, 0, 1), clamp(c.g, 0, 1), clamp(c.b, 0, 1)) for c in permutedims(img[:, end:-1:1])]

"Halve an RGB image by 2×2 averaging (the last row/column is repeated for odd sizes)."
function halve(a)
    nr, nc = size(a)
    out = Matrix{RGB{N0f8}}(undef, cld(nr, 2), cld(nc, 2))
    for j in axes(out, 2), i in axes(out, 1)
        r0, c0 = 2i - 1, 2j - 1
        ps = (a[r0, c0], a[min(r0 + 1, nr), c0], a[r0, min(c0 + 1, nc)], a[min(r0 + 1, nr), min(c0 + 1, nc)])
        out[i, j] = RGB{N0f8}(sum(float(p.r) for p in ps)/4, sum(float(p.g) for p in ps)/4, sum(float(p.b) for p in ps)/4)
    end
    return out
end

"Number of zoom levels below the native one such that level 0 fits in about one tile."
zmax_for(n) = max(0, ceil(Int, log2(n/TILE)))

"Write the tile pyramid of an image (rows top-down) to dir/<z>/<x>/<y>.jpg; returns zmax."
function write_tiles(a, dir; quality=85)
    zmax = zmax_for(maximum(size(a)))
    for z in zmax:-1:0
        nr, nc = size(a)
        for tx in 0:cld(nc, TILE)-1, ty in 0:cld(nr, TILE)-1
            t = fill(RGB{N0f8}(1, 1, 1), TILE, TILE)            # white beyond the grid edge
            rr = ty*TILE+1:min((ty + 1)*TILE, nr); cc = tx*TILE+1:min((tx + 1)*TILE, nc)
            t[1:length(rr), 1:length(cc)] .= a[rr, cc]
            f = joinpath(dir, string(z), string(tx), "$(ty).jpg")
            mkpath(dirname(f))
            jpeg_encode(f, t; quality)
        end
        z > 0 && (a = halve(a))
    end
    return zmax
end

# ---------------------------------------------------------------------------
# GeoJSON (projected coordinates in metres)
# ---------------------------------------------------------------------------

"Split NaN-separated vertices (km) into line strings (m), dropping single points."
function linestrings(pts)
    out = Vector{Vector{Vector{Float64}}}(); cur = Vector{Vector{Float64}}()
    for p in vcat(pts, [Point2f(NaN, NaN)])
        if isnan(p[1])
            length(cur) >= 2 && push!(out, cur); cur = Vector{Vector{Float64}}()
        else
            push!(cur, [round(1e3p[1]; digits=-1), round(1e3p[2]; digits=-1)])
        end
    end
    return out
end

line_feature(pts, props) = Dict("type" => "Feature", "properties" => props,
                                "geometry" => Dict("type" => "MultiLineString", "coordinates" => linestrings(pts)))
point_feature(xy, props) = Dict("type" => "Feature", "properties" => props,
                                "geometry" => Dict("type" => "Point", "coordinates" => [round(1e3xy[1]), round(1e3xy[2])]))

function write_geojson(path, features)
    mkpath(dirname(path))
    open(io -> JSON.json(io, Dict("type" => "FeatureCollection", "features" => features)), path, "w")
end

# full names of the IMBIE 2 Greenland regions, searchable on the web map
const REGION_FULL = Dict("NO" => "North", "NE" => "Northeast", "SE" => "Southeast", "SW" => "Southwest",
                         "CW" => "Central west", "NW" => "Northwest")

label_props(l::PlaceLabel) = Dict("name" => l.name, "type" => l.type, "tier" => l.tier, "source" => l.source,
                                  "alt" => join(l.alt, ", "), "lat" => round(l.lat, digits=3), "lon" => round(l.lon, digits=3))

# ---------------------------------------------------------------------------
# Colour bars for the legend (CSS gradients)
# ---------------------------------------------------------------------------

hexcolor(c) = "#" * hex(RGB{N0f8}(c))
gradient_stops(cs; n=16) = [hexcolor(get(cs, t)) for t in range(0, 1, length=n)]

function colorbars(style, velcmap, lim)
    vel = Dict("label" => "Surface ice velocity [m/yr]", "stops" => gradient_stops(VEL_CMAPS[velcmap]),
               "ticks" => [[(log10(v) - 0.3)/3.2, string(v)] for v in (10, 100, 1000)])
    srf = Dict("label" => "Ice-surface elevation [m]", "stops" => gradient_stops(CS_ICE),
               "ticks" => [[v/lim.srflim[2], string(v)] for v in 0:1000:lim.srflim[2]])
    bed = Dict("label" => "Bed elevation [m]",
               "stops" => [hexcolor(bedcolor(v, lim.bedlim...)) for v in range(lim.bedlim..., length=24)],
               "ticks" => [[(v - lim.bedlim[1])/(lim.bedlim[2] - lim.bedlim[1]), string(v)] for v in -1000:1000:3000
                           if lim.bedlim[1] <= v <= lim.bedlim[2]])
    ocean = Dict("label" => "Ocean depth [m]", "stops" => gradient_stops(reverse(CS_OCEAN_LIGHT)),
                 "ticks" => [[v/4500, string(v)] for v in 0:1000:4000])
    style === :velocity && return [vel, ocean]
    style === :surface  && return [srf, ocean]
    return [bed]
end

# ---------------------------------------------------------------------------

function build_region(region; full)
    out = joinpath(ASSETS, region); mkpath(out)
    proj = PROJS[region]; g = GRIDS[region]; lim = STYLE_LIMITS[region]
    d = load_prepared(region; full)
    r = d.r; x, y, dx = r.x, r.y, r.dx

    # raster tiles of each style
    zmax = 0
    for (id, _, style, cm) in WEB_STYLES
        img = compose(style, r; lim..., ocean=:light, velcmap=VEL_CMAPS[cm])
        dir = joinpath(out, "tiles", id); rm(dir; recursive=true, force=true)
        zmax = write_tiles(rows_top_down(img), dir)
        println("tiles: ", region, " ", id, " (", zmax + 1, " levels)")
    end

    # vector layers; contours from a grid of ≥ 1 km spacing to keep the files small
    st = max(1, round(Int, 1.0/dx))
    xs, ys = x[1:st:end], y[1:st:end]
    levels = 500:500:(region == "greenland" ? 3000 : 4000)
    write_geojson(joinpath(out, "contours.geojson"),
        [line_feature(elevation_contour_points(xs, ys, r.zs[1:st:end, 1:st:end], d.grounded[1:st:end, 1:st:end], [lv]),
                      Dict("elevation" => lv)) for lv in levels])
    write_geojson(joinpath(out, "margin.geojson"),
        [line_feature(maskoutline_points(xs, ys, r.ice[1:st:end, 1:st:end]; σ=1.0), Dict("kind" => "ice margin")),
         line_feature(maskoutline_points(xs, ys, d.grounded[1:st:end, 1:st:end]; σ=1.0, mask=r.ice[1:st:end, 1:st:end]),
                      Dict("kind" => "grounding line"))])
    icemask_b = region == "greenland" ? d.ice_b : d.grounded_b
    write_geojson(joinpath(out, "divides.geojson"),
        [line_feature(basin_divide_points(d.xb, d.yb, d.basin, icemask_b; σ=1.2), Dict())])
    write_geojson(joinpath(out, "regions.geojson"),
        [point_feature(b.p, Dict("name" => b.name, "type" => "drainage region", "alt" => get(REGION_FULL, b.name, ""),
                                 "extent" => [round(1e3v) for v in b.extent]))
         for b in basin_label_points(d.xb, d.yb, d.basin, d.basin_names, d.ice_b)])
    labs = read_labels(joinpath(ROOT, "data", "labels_$(region).csv"))
    gaz  = read_labels(joinpath(REPO_PREP_DIR, "gazetteer_$(region).csv"))
    write_geojson(joinpath(out, "labels.geojson"), [point_feature(project(proj, l.lon, l.lat), label_props(l)) for l in labs])
    write_geojson(joinpath(out, "gazetteer.geojson"), [point_feature(project(proj, l.lon, l.lat), label_props(l)) for l in gaz])
    seaice = []
    for (mm, dash) in zip(SEAICE_MONTHS[region], ([10, 8], [2, 6]))
        f = "seaice_$(mm).geojson"
        fc = JSON.parsefile(joinpath(REPO_PREP_DIR, "seaice_$(region)_$(mm).geojson"))
        write_geojson(joinpath(out, f), fc["features"])        # without the crs member (set by the map page)
        push!(seaice, Dict("file" => f, "dash" => dash, "label" => "Sea-ice edge, $(MONTHS[parse(Int, mm)]) (median 1981–2010)"))
    end

    # configuration of the map page
    x0, x1 = 1e3*(x[1] - dx/2), 1e3*(x[end] + dx/2)
    y0, y1 = 1e3*(y[1] - dx/2), 1e3*(y[end] + dx/2)
    a = d.attrs
    cfg = Dict(
        "region" => region, "epsg" => "EPSG:$(g.epsg)", "proj4" => PROJ4[region],
        "extent" => [x0, y0, x1, y1], "origin" => [x0, y1], "tileSize" => TILE,
        "resolutions" => [1e3*dx*2.0^(zmax - z) for z in 0:zmax],
        "styles" => [Dict("id" => id, "name" => name, "colorbars" => colorbars(style, cm, lim))
                     for (id, name, style, cm) in WEB_STYLES],
        "seaice" => seaice,
        "numbers" => Dict(String(k) => Float64(v) for (k, v) in a if v isa Real),
        "resolution_m" => 1e3dx)
    open(io -> JSON.json(io, cfg), joinpath(out, "config.json"), "w")
    println("wrote ", out)
end

# readable descriptions of the parts of poster file names
const NAME_PARTS = Dict("velocity" => "ice velocity", "surface" => "surface elevation", "bed" => "bed topography",
    "dark" => "dark ocean", "nocontours" => "no contours", "ember" => "ember colours", "batlow" => "batlow colours",
    "lajolla" => "lajolla colours", "tier2" => "more names", "custom" => "added names",
    "vintage" => "old-world style")

const REGION_NAMES = Dict("greenland" => "Greenland", "antarctica" => "Antarctica", "nh" => "Northern Hemisphere")

"Caption from a poster file name: region[, time slice]: options, e.g. \"Northern Hemisphere, 20 ka: surface elevation\"."
function caption(f)
    parts = filter(!=("A0"), split(replace(f, "_small.png" => ""), "_"))
    name = REGION_NAMES[popfirst!(parts)]
    occursin(r"^[0-9.]+ka$", first(parts)) && (name *= ", " * replace(popfirst!(parts), "ka" => " ka"))
    return name * ": " * join([get(NAME_PARTS, p, p) for p in parts], ", ")
end

"""
Gallery section `title` of the small PNGs in `dir` (`cols` of 12 grid columns
each), copied to assets/img/, with links to the PDF and full PNG where they
exist (copied to assets/posters/), or else to the small PNG itself (`pnglink`).
"""
function gallery_section!(md, dir, title; cols=4, pnglink=false)
    img, big = joinpath(ASSETS, "img"), joinpath(ASSETS, "posters")
    isdir(dir) || return
    fs = sort(filter(f -> endswith(f, "_small.png") && !occursin("_coarse", f), readdir(dir)); rev=true)
    isempty(fs) && return
    println(md, "## ", title, "\n\n::: {.grid .poster-grid}")
    for f in fs
        cp(joinpath(dir, f), joinpath(img, f); force=true)
        println(md, "::: {.g-col-12 .g-col-md-", cols, "}")
        println(md, "![", caption(f), "](assets/img/", f, "){group=\"", title, "\"}\n")
        links = String[]
        for (ext, what) in ((".pdf", "PDF, A0"), (".png", "PNG, 100 dpi"))
            src = joinpath(dir, replace(f, "_small.png" => ext))
            isfile(src) || continue
            cp(src, joinpath(big, basename(src)); force=true)
            mb = round(Int, filesize(src)/2^20)
            push!(links, "[$what, $mb MB](assets/posters/$(basename(src)))")
        end
        if isempty(links) && pnglink
            push!(links, "[PNG, 1600 px, $(round(Int, filesize(joinpath(dir, f))/2^20)) MB](assets/img/$f)")
        end
        println(md, "[", caption(f), "]{.poster-links}", isempty(links) ? "" : " · " * join(links, " · "))
        println(md, ":::")
    end
    println(md, ":::\n")
end

"""
Copy the poster files to the site (small PNGs to assets/img/, the PDFs and
full PNGs of the default posters to assets/posters/) and write the galleries:
site/_gallery.md (included by index.qmd) with the default posters and their
downloads, then the variants, site/_gallery_paleo.md (paleo.qmd) with the
paleo maps and site/_gallery_vintage.md (historical.qmd) with the vintage maps.
"""
function build_gallery()
    for d in (joinpath(ASSETS, "img"), joinpath(ASSETS, "posters"))
        rm(d; recursive=true, force=true); mkpath(d)
    end
    md = IOBuffer()
    gallery_section!(md, joinpath(ROOT, "plots"), "Posters"; cols=6)
    gallery_section!(md, joinpath(ROOT, "plots", "variants"), "Variants")
    write(joinpath(SITE, "_gallery.md"), take!(md))
    gallery_section!(md, joinpath(ROOT, "plots", "paleo"), "Maps"; cols=6, pnglink=true)
    write(joinpath(SITE, "_gallery_paleo.md"), take!(md))
    gallery_section!(md, joinpath(ROOT, "plots", "vintage"), "Maps"; cols=6)
    write(joinpath(SITE, "_gallery_vintage.md"), take!(md))
end

if abspath(PROGRAM_FILE) == @__FILE__
    regions = filter(in(("greenland", "antarctica")), ARGS)
    if !("gallery" in ARGS)
        for reg in (isempty(regions) ? ["greenland", "antarctica"] : regions)
            build_region(reg; full="full" in ARGS)
        end
    end
    build_gallery()
end
