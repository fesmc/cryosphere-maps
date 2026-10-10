# Old-world style A0 poster maps: the same data as the posters of greenland.jl and
# antarctica.jl (ice velocity, topography, sea-ice edge, place names) drawn in
# inks on paper, with water lining along the coasts, waves on the open sea, pack
# ice inside the winter sea-ice edge, expedition routes and wildlife.
# Usage: julia --project=. scripts/vintage.jl antarctica|greenland [noroutes] [nofauna] [full] [pdf]
# `full` and `pdf` as for the other posters (see parse_args in common.jl).

include("common.jl")
using Random

# fonts: IM Fell English (Igino Marini, SIL Open Font License, data/fonts/OFL.txt)
const FONT_DIR = joinpath(ROOT, "data", "fonts")
const VF = (regular=joinpath(FONT_DIR, "IMFeENrm28P.ttf"), italic=joinpath(FONT_DIR, "IMFeENit28P.ttf"),
            caps=joinpath(FONT_DIR, "IMFeENsc28P.ttf"))

# paper and inks
const PAPER    = RGBf(colorant"#efe2c2")
const V_SEA    = RGBf(colorant"#cfd3bd")
const V_DEEP   = RGBf(colorant"#9fb0a6")
const V_PACK   = RGBf(colorant"#e9e8dc")
const V_ICE    = RGBf(colorant"#f7f0dc")
const V_SHELF  = RGBf(colorant"#e6e3d0")
const V_ROCK   = RGBf(colorant"#8a6a48")
const INK      = colorant"#3b2a1a"
const SEPIA    = colorant"#6b4a2b"
const SEA_INK  = colorant"#2f4f5f"
const RED_INK  = colorant"#8e1b1b"
const CS_FLOW  = cgrad([colorant"#ecd9a8", colorant"#d1a462", colorant"#b06a3b", colorant"#8e3a2a", colorant"#6a2330"])
const FLOW_LIM = (0.6, 3.3)                    # log10 speed range of CS_FLOW
const V_ROUTE_COLORS = [RED_INK, colorant"#1f3f6b", colorant"#4d5d24", colorant"#6a2c5e", colorant"#5a3a1a",
                        colorant"#2e6b6b", colorant"#a0522d", colorant"#4b3b7a"]

"Sea names in spaced capitals, one word per line for names of more than `maxlen` letters."
function spaced(s; maxlen=12)
    words = [join(collect(uppercase(w)), ' ') for w in split(s)]
    return join(words, count(isletter, s) > maxlen ? "\n" : "   ")
end

const LSTYLE_VINTAGE = Dict(
    "region"     => (size=22, font=VF.caps,    color=INK,     case=uppercase, outer=false, marker=false),
    "mountains"  => (size=17, font=VF.italic,  color=SEPIA,   case=identity,  outer=false, marker=false),
    "basin"      => (size=15, font=VF.italic,  color=SEPIA,   case=identity,  outer=false, marker=false),
    "dome"       => (size=14, font=VF.regular, color=INK,     case=identity,  outer=false, marker=true),
    "station"    => (size=14, font=VF.regular, color=INK,     case=identity,  outer=false, marker=true),
    "icecore"    => (size=14, font=VF.regular, color=RED_INK, case=identity,  outer=false, marker=true),
    "lake"       => (size=14, font=VF.italic,  color=SEA_INK, case=identity,  outer=false, marker=true),
    "island"     => (size=14, font=VF.regular, color=INK,     case=identity,  outer=false, marker=false),
    "sea"        => (size=22, font=VF.italic,  color=SEA_INK, case=spaced,    outer=false, marker=false),
    "settlement" => (size=14, font=VF.regular, color=INK,     case=identity,  outer=true,  marker=true),
    "glacier"    => (size=14, font=VF.italic,  color=INK,     case=identity,  outer=true,  marker=false),
    "iceshelf"   => (size=14, font=VF.italic,  color=INK,     case=identity,  outer=true,  marker=false),
    "fjord"      => (size=14, font=VF.italic,  color=SEA_INK, case=identity,  outer=true,  marker=false),
    "site"       => (size=13, font=VF.italic,  color=INK,     case=identity,  outer=false, marker=true),
    "pole"       => (size=16, font=VF.italic,  color=RED_INK, case=identity,  outer=false, marker=true),
)

# ---------------------------------------------------------------------------
# Raster and sea
# ---------------------------------------------------------------------------

"Smooth multi-scale noise in about [-1, 1], for the mottling of the paper."
function paper_noise(nx, ny; seed=3)
    rng = MersenneTwister(seed)
    n = zeros(nx, ny)
    for (σ, w) in ((40.0, 0.5), (10.0, 0.3), (2.0, 0.2))
        a = gauss_smooth(randn(rng, nx, ny), σ)
        n .+= w .* a ./ maximum(abs, a)
    end
    return n
end

"""
RGB image in inks on paper: the sea tinted by depth (lighter where `pack` marks
the winter pack ice), ice in cream with the hillshade, fast ice in rust to
oxblood, rock in sepia.
"""
function compose_vintage(r, pack)
    nx, ny = size(r.zs)
    nz = paper_noise(nx, ny)
    img = Matrix{RGBf}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        ice, ocean = r.ice[i, j], r.ocean[i, j] && !r.ice[i, j]
        if ocean
            c = mix(mix(PAPER, V_SEA, 0.55), V_DEEP, 0.55*clamp(-r.zb[i, j]/4500, 0, 1))
            pack[i, j] && (c = mix(c, V_PACK, 0.45))
            c = shade(c, r.hs_b[i, j]; strength=0.25)
        elseif !ice
            c = shade(V_ROCK, r.hs_b[i, j]; strength=1.2)
        else
            c = shade(r.shelf[i, j] ? V_SHELF : V_ICE, r.hs_s[i, j]; strength=0.9)
            u = r.u[i, j]
            if !isnan(u) && u > 0
                lu = log10(max(u, 1.0))
                c = mix(c, cmap(CS_FLOW, lu, FLOW_LIM...), 0.75*clamp((lu - 1.0)/1.2, 0, 1))
            end
        end
        img[i, j] = shade(mix(c, PAPER, 0.12), 0.06nz[i, j])
    end
    return img
end

"Chamfer distance (km) from the `src` cells on a grid of spacing dx km."
function distance_km(src, dx)
    nx, ny = size(src)
    d = [s ? 0.0 : Inf for s in src]
    a, b = 1.0, sqrt(2)
    for j in 1:ny, i in 1:nx
        i > 1 && (d[i, j] = min(d[i, j], d[i-1, j] + a))
        j > 1 && (d[i, j] = min(d[i, j], d[i, j-1] + a))
        i > 1 && j > 1 && (d[i, j] = min(d[i, j], d[i-1, j-1] + b))
        i < nx && j > 1 && (d[i, j] = min(d[i, j], d[i+1, j-1] + b))
    end
    for j in ny:-1:1, i in nx:-1:1
        i < nx && (d[i, j] = min(d[i, j], d[i+1, j] + a))
        j < ny && (d[i, j] = min(d[i, j], d[i, j+1] + a))
        i < nx && j < ny && (d[i, j] = min(d[i, j], d[i+1, j+1] + b))
        i > 1 && j < ny && (d[i, j] = min(d[i, j], d[i-1, j+1] + b))
    end
    return d .* dx
end

gridindex(r, p) = (clamp(round(Int, (p[1] - r.x[1])/r.dx) + 1, 1, length(r.x)),
                   clamp(round(Int, (p[2] - r.y[1])/r.dx) + 1, 1, length(r.y)))

"Split NaN-separated line vertices into their parts."
function line_parts(pts)
    parts = [Point2f[]]
    for p in pts
        isnan(p[1]) ? (isempty(parts[end]) || push!(parts, Point2f[])) : push!(parts[end], p)
    end
    return filter(!isempty, parts)
end

"""
Winter pack ice: open-water cells (ocean, no ice shelf) that cannot be reached
from the `seeds` (km, points on the open sea) without crossing the median
sea-ice edge `edge` (NaN-separated line vertices, km). The edge is split where
it crosses land; ends of parts that stop inside the grid are joined to the
nearest end of another part within `maxgap` km, so the edge stays closed
across straits narrower than the grid resolves.
"""
function pack_ice(r, edge, seeds; maxgap=300.0)
    nx, ny = size(r.zs)
    water = r.ocean .& .!r.ice
    wall = falses(nx, ny)
    draw(a, b) = for t in range(0, 1, length=ceil(Int, 2hypot((b - a)...)/r.dx) + 2)
        wall[gridindex(r, a + t*(b - a))...] = true
    end
    parts = line_parts(edge)
    for p in parts, k in 1:length(p)-1
        draw(p[k], p[k+1])
    end
    m = 2r.dx
    inner(q) = r.x[1] + m < q[1] < r.x[end] - m && r.y[1] + m < q[2] < r.y[end] - m
    ends = [(k, q) for (k, p) in enumerate(parts) for q in (p[1], p[end]) if inner(q)]
    for (k, q) in ends
        others = [e for (j, e) in ends if j != k]
        isempty(others) && continue
        e = others[argmin([hypot((e - q)...) for e in others])]
        hypot((e - q)...) <= maxgap && draw(q, e)
    end
    reach = falses(nx, ny)
    stack = [gridindex(r, s) for s in seeds]
    while !isempty(stack)
        i, j = pop!(stack)
        (1 <= i <= nx && 1 <= j <= ny) || continue
        (reach[i, j] || wall[i, j] || !water[i, j]) && continue
        reach[i, j] = true
        push!(stack, (i+1, j), (i-1, j), (i, j+1), (i, j-1))
    end
    return water .& .!reach
end

const G_WAVE = BezierPath("M0,5 C2,0 4,0 6,5 C8,10 10,10 12,5 C14,0 16,0 18,5"; fit=true, flipy=true)

"Ink lines offsetting the coast at `dists` (km), fading seaward."
function water_lines!(ax, xs, ys, dist; dists=(10, 22, 36, 52, 72, 96))
    for (k, dkm) in enumerate(dists)
        lines!(ax, contour_points(xs, ys, dist, [dkm]); color=(INK, 0.75 - 0.1k), linewidth=1.6 - 0.18k)
    end
end

"""
Wave strokes on a jittered grid over open water beyond the water lining, and
floes scattered over the pack ice. Returns nothing.
"""
function sea_marks!(ax, r, dist, pack, limits; spacing=110.0, nfloes=4000, wavesize=26, seed=7)
    rng = MersenneTwister(seed)
    open = r.ocean .& .!r.ice .& .!pack
    waves = Point2f[]
    for (row, y0) in enumerate(limits[3]:0.8spacing:limits[4]), x0 in limits[1]:spacing:limits[2]
        p = (x0 + (isodd(row) ? spacing/2 : 0) + 0.25spacing*randn(rng), y0 + 0.2spacing*randn(rng))
        i, j = gridindex(r, p)
        open[i, j] && dist[i, j] > 110 && push!(waves, Point2f(p))
    end
    scatter!(ax, waves; marker=G_WAVE, markersize=wavesize, color=(SEPIA, 0.55))
    floes = Point2f[]
    for _ in 1:nfloes
        p = (limits[1] + rand(rng)*(limits[2] - limits[1]), limits[3] + rand(rng)*(limits[4] - limits[3]))
        i, j = gridindex(r, p)
        pack[i, j] && dist[i, j] > 30 && push!(floes, Point2f(p))
    end
    scatter!(ax, floes; marker=:hexagon, markersize=10, color=(:white, 0.9), strokecolor=(INK, 0.55), strokewidth=0.9)
end

# ---------------------------------------------------------------------------
# Map furniture
# ---------------------------------------------------------------------------

"Map angle (degrees clockwise from up) of true north at the projected point (x, y) km."
function north_angle(proj, x, y)
    lon, lat = unproject(proj, x, y)
    q = project(proj, lon, lat + 0.1)
    return rad2deg(atan(q[1] - x, q[2] - y))
end

"Eight-point compass rose centred at (x, y) km with radius R km; `north` is the map angle of north (degrees)."
function compass_rose!(ax, x, y, R; north=0.0, fontsize=40)
    for (k, θ) in enumerate(0:45:315)
        L = iseven(k) ? 0.55R : R; w = iseven(k) ? 0.09R : 0.14R
        φ = deg2rad(θ + north)
        tip = Point2f(x + L*sin(φ), y + L*cos(φ))
        l  = Point2f(x + w*sin(φ - π/2), y + w*cos(φ - π/2))
        rt = Point2f(x + w*sin(φ + π/2), y + w*cos(φ + π/2))
        poly!(ax, [Point2f(x, y), l, tip]; color=INK, strokecolor=INK, strokewidth=0.8)
        poly!(ax, [Point2f(x, y), tip, rt]; color=PAPER, strokecolor=INK, strokewidth=0.8)
    end
    for rr in (0.62R, 0.68R)
        lines!(ax, [Point2f(x + rr*cos(t), y + rr*sin(t)) for t in range(0, 2π, length=200)]; color=INK, linewidth=1.2)
    end
    φ = deg2rad(north)
    text!(ax, x + 1.18R*sin(φ), y + 1.18R*cos(φ); text="N", font=VF.caps, fontsize, color=INK,
          align=(:center, :center), rotation=-φ)
    return (x - 1.3R, x + 1.3R, y - 1.3R, y + 1.3R)
end

"Double neatline with alternating ink bars every `bar` km, width `b` km. Returns the frame obstacles."
function neatline!(ax, limits; bar=200.0, b=25.0)
    x0, x1, y0, y1 = limits
    lines!(ax, Rect(x0, y0, x1 - x0, y1 - y0); color=INK, linewidth=3)
    lines!(ax, Rect(x0 + b, y0 + b, x1 - x0 - 2b, y1 - y0 - 2b); color=INK, linewidth=1.2)
    for (k, s) in enumerate(x0:bar:x1 - bar), yy in (y0, y1 - b)
        isodd(k) && poly!(ax, Rect(s, yy, bar, b); color=INK)
    end
    for (k, s) in enumerate(y0:bar:y1 - bar), xx in (x0, x1 - b)
        isodd(k) && poly!(ax, Rect(xx, s, b, bar); color=INK)
    end
    return [(x0, x1, y0, y0 + b), (x0, x1, y1 - b, y1), (x0, x0 + b, y0, y1), (x1 - b, x1, y0, y1)]
end

"Ice-flow scale in the vintage inks."
flow_colorbar!(pos; width, kw...) =
    Colorbar(pos; colormap=CS_FLOW, limits=FLOW_LIM, vertical=false, width, height=26, flipaxis=true,
             ticks=([1, 2, 3], ["10", "100", "1000"]), ticklabelfont=VF.regular, ticklabelsize=26,
             ticklabelcolor=INK, label="Ice flow, metres a year", labelfont=VF.italic, labelsize=28,
             labelcolor=INK, spinewidth=1, topspinecolor=INK, bottomspinecolor=INK, leftspinecolor=INK,
             rightspinecolor=INK, kw...)

vintage_legend!(pos, entries; kw...) =
    Legend(pos, first.(entries), last.(entries); framevisible=false, backgroundcolor=:transparent,
           labelfont=VF.regular, labelsize=26, labelcolor=INK, patchsize=(50, 22), rowgap=4, kw...)

# ---------------------------------------------------------------------------
# The map
# ---------------------------------------------------------------------------

"""
Everything inside the map axis: raster, water lining, sea marks, outlines,
graticule, routes, labels, wildlife, frame, compass rose and scale bar.
Returns the legend entries.
"""
function vintage_map!(ax, d, region; limits, kmpp, scale, proj, labels, routes, fauna, seaice, seeds,
                      graticule, rose, scalebar, bar, faunanames=Dict{String, String}())
    r = d.r
    ix = findall(xi -> limits[1] <= xi <= limits[2], r.x); iy = findall(yi -> limits[3] <= yi <= limits[4], r.y)
    xs = r.x[ix]; ys = r.y[iy]; h = r.dx/2
    edge = geojson_lines(joinpath(REPO_PREP_DIR, "seaice_$(region)_$(seaice).geojson"))
    pack = pack_ice(r, edge, seeds)
    img = compose_vintage(r, pack)
    image!(ax, (xs[1] - h, xs[end] + h), (ys[1] - h, ys[end] + h), img[ix, iy]; interpolate=true)

    dist = gauss_smooth(distance_km(.!(r.ocean .& .!r.ice), r.dx), 1.5)
    water_lines!(ax, xs, ys, dist[ix, iy])
    sea_marks!(ax, r, dist, pack, limits)
    lines!(ax, edge; color=(SEA_INK, 0.7), linewidth=2.0, linestyle=:dash)

    maskoutline!(ax, xs, ys, r.ice[ix, iy]; σ=1.5, color=INK, linewidth=2.2)
    maskoutline!(ax, xs, ys, d.grounded[ix, iy]; σ=1.5, mask=r.ice[ix, iy], color=(INK, 0.6), linewidth=1.2)
    elevation_contours!(ax, xs, ys, r.zs[ix, iy], d.grounded[ix, iy], 1000:1000:4000; color=(SEPIA, 0.5), linewidth=1.0)
    inside = p -> limits[1] <= p[1] <= limits[2] && limits[3] <= p[2] <= limits[4]
    graticule!(ax, proj, graticule.lats, graticule.lons; graticule.latrange, color=(SEPIA, 0.45), lw=1.0, inside)
    circle = [Point2f(project(proj, λ, sign(proj.lat_ts)*66.56)...) for λ in -180:0.5:180]
    lines!(ax, [inside(p) ? p : Point2f(NaN, NaN) for p in circle]; color=(RED_INK, 0.6), linewidth=1.6)

    frame = neatline!(ax, limits; bar)
    obstacles = vcat(frame, [compass_rose!(ax, rose...; north=north_angle(proj, rose[1], rose[2]))])
    scalebar!(ax, scalebar..., region == "greenland" ? 400.0 : 1000.0; fontsize=26, color=INK, font=VF.regular,
              ink=INK, paper=PAPER)
    url = site_url!(ax, limits, kmpp; corner=:right, fontsize=24, margin=40, font=VF.italic, color=(INK, 0.7))
    push!(obstacles, url)

    rbox, sites, rlegend = routes!(ax, routes, proj; kmpp, colors=V_ROUTE_COLORS, scale, font=VF.italic, paper=PAPER,
                                   obstacles, limits)
    lbox = draw_labels!(ax, vcat(labels, sites), proj; layout=:coastal, kmpp, scale, styles=LSTYLE_VINTAGE,
                        ink=INK, paper=PAPER, seacolor=SEA_INK, obstacles=vcat(obstacles, rbox), limits,
                        icemask=r.ice, shelfmask=r.shelf, groundedmask=d.grounded, x=r.x, y=r.y,
                        offset=region == "greenland" ? 50.0 : 120.0, maxlead=region == "greenland" ? 300.0 : 500.0,
                        tmax=region == "greenland" ? 120.0 : 250.0)
    _, flegend = fauna!(ax, fauna, proj; kmpp, scale, ink=INK, paper=V_ICE, obstacles=vcat(obstacles, rbox, lbox),
                        names=faunanames)

    return vcat([(LineElement(color=INK, linewidth=2.2), "Ice front"),
                 (LineElement(color=(INK, 0.6), linewidth=1.2), "Grounding line"),
                 (LineElement(color=(SEPIA, 0.5), linewidth=1.0), "1000 m surface contours"),
                 (LineElement(color=(SEA_INK, 0.7), linewidth=2.0, linestyle=:dash),
                  "Sea-ice edge, $(MONTHS[parse(Int, seaice)])"),
                 (MarkerElement(marker=G_WAVE, color=(SEPIA, 0.7), markersize=30), "Open sea in winter"),
                 (MarkerElement(marker=:hexagon, color=(:white, 0.9), strokecolor=(INK, 0.55), strokewidth=0.9,
                                markersize=12), "Pack ice in winter")],
                rlegend, flegend)
end

const CREDITS_VINTAGE = "Routes from the expedition accounts and published positions (data/routes_REGION.csv). "

# Greenlandic (kalaallisut) words on the Greenland map
const GREENLANDIC_FAUNA = Dict("polarbear" => "Nanoq · polar bear denning area",
                               "muskox"    => "Umimmak · musk-ox calving area",
                               "walrus"    => "Aaveq · walrus haul-out",
                               "narwhal"   => "Qilalugaq qernertaq · narwhal summer area")
const GREENLANDIC_GLOSSARY = "sermersuaq the ice sheet · sermeq glacier · kangerluk fjord\n" *
                             "qeqertaq island · imaq sea · nuna land"

# ---------------------------------------------------------------------------
# Posters
# ---------------------------------------------------------------------------

function plot_vintage_antarctica(d; labels, routes, fauna, scale=1.85)
    xmap = (-3040.0, 3040.0); ymap = (-2600.0, 2500.0); limits = (xmap..., ymap...)
    pad = 70
    H = A0_LANDSCAPE[2] - 2pad; kmpp = (ymap[2] - ymap[1])/H; W = (xmap[2] - xmap[1])/kmpp
    panel = A0_LANDSCAPE[1] - 2pad - W - 60

    fig = Figure(size=A0_LANDSCAPE, figure_padding=pad, backgroundcolor=PAPER, fontsize=24)
    ax = Axis(fig[1, 1]; width=W, height=H, limits=(xmap, ymap), backgroundcolor=PAPER)
    hidedecorations!(ax); hidespines!(ax)
    # fewer names than on the velocity poster: no glaciers, ice streams, basins or lakes
    keep(l) = l.type ∉ ("glacier", "basin", "lake") && l.name != "Amundsen-Scott South Pole Station"
    labs = vcat(filter(keep, labels),
                PlaceLabel("Polus Australis", 0.0, -90.0, "pole", 1, "", String[]))
    legend = vintage_map!(ax, d, "antarctica"; limits, kmpp, scale, proj=PROJ_ANT, labels=labs, routes, fauna,
                          seaice="09", seeds=[(-3000.0, 3000.0), (3000.0, 3000.0)],
                          graticule=(lats=-85:5:-60, lons=-180:15:165, latrange=(-88, -55)),
                          rose=(-2420.0, -1780.0, 260.0), scalebar=(xmap[1] + 200, ymap[1] + 180), bar=200.0)

    side = fig[1, 2] = GridLayout(width=panel)
    Label(side[1, 1], "Terra\nAustralis"; font=VF.caps, fontsize=110, color=INK, lineheight=0.9)
    Label(side[2, 1], "The Antarctic Ice Sheet"; font=VF.italic, fontsize=46, color=SEPIA)
    Label(side[3, 1], "its ice in motion, the voyages of\ndiscovery & the rookeries of\nthe penguins";
          font=VF.italic, fontsize=30, color=INK, justification=:center)
    flow_colorbar!(side[4, 1]; width=panel - 60)
    vintage_legend!(side[5, 1], legend; halign=:left)
    qr_code!(side[6, 1]; size=130, fontsize=20, ink=INK, paper=PAPER, font=VF.italic, color=INK)
    Label(side[7, 1], credits_vintage("antarctica"); font=VF.regular, fontsize=19, color=(INK, 0.8),
          word_wrap=true, width=panel, justification=:left)
    rowgap!(side, 1, 10); rowgap!(side, 3, 50); rowgap!(side, 4, 40); rowgap!(side, 5, 40)
    colgap!(fig.layout, 60)
    return fig
end

function plot_vintage_greenland(d; labels, routes, fauna, scale=1.85)
    pad = 70
    W = A0_PORTRAIT[1] - 2pad; H = 2560.0
    ymap = (-3420.0, -600.0); kmpp = (ymap[2] - ymap[1])/H
    xc = 120.0; xl = (xc - W*kmpp/2, xc + W*kmpp/2); limits = (xl..., ymap...)

    fig = Figure(size=A0_PORTRAIT, figure_padding=pad, backgroundcolor=PAPER, fontsize=24)
    head = fig[1, 1] = GridLayout()
    Label(head[1, 1], "Kalaallit Nunaat"; font=VF.caps, fontsize=120, color=INK)
    Label(head[2, 1], "Sermersuaq · the Greenland Ice Sheet: its ice in motion, the routes of its crossings and its wildlife";
          font=VF.italic, fontsize=40, color=SEPIA)
    rowgap!(head, 0)
    ax = Axis(fig[2, 1]; width=W, height=H, limits=(xl, ymap), backgroundcolor=PAPER)
    hidedecorations!(ax); hidespines!(ax)
    legend = vintage_map!(ax, d, "greenland"; limits, kmpp, scale, proj=PROJ_GRL, labels, routes, fauna,
                          seaice="03", seeds=[(limits[2] - 50, limits[3] + 50)],
                          graticule=(lats=60:5:80, lons=-80:10:0, latrange=(58, 84)),
                          rose=(xl[1] + 260.0, ymap[1] + 420.0, 150.0), scalebar=(700.0, ymap[1] + 110), bar=100.0,
                          faunanames=GREENLANDIC_FAUNA)

    foot = fig[3, 1] = GridLayout()
    flow_colorbar!(foot[1, 1]; width=500, valign=:top)
    Label(foot[2, 1], GREENLANDIC_GLOSSARY; font=VF.italic, fontsize=22, color=INK, justification=:left,
          halign=:left, valign=:top, word_wrap=true, width=500)
    vintage_legend!(foot[1, 2], legend; nbanks=3, valign=:top)
    qr_code!(foot[1, 3]; size=130, fontsize=20, ink=INK, paper=PAPER, font=VF.italic, color=INK, valign=:top)
    colgap!(foot, 50)
    Label(fig[4, 1], credits_vintage("greenland"); font=VF.regular, fontsize=21, color=(INK, 0.8), word_wrap=true,
          width=W - 200, justification=:left)
    rowgap!(fig.layout, 25)
    return fig
end

"Data credits of the vintage maps."
function credits_vintage(region)
    base = region == "antarctica" ?
        "Topography from BedMachine Antarctica v4 (Morlighem et al., 2020); ice velocity from MEaSUREs InSAR v2 " *
        "(Rignot et al., 2011); median sea-ice edge 1981–2010 from the Sea Ice Index v4 (Fetterer et al., 2017); " *
        "place names from the SCAR Composite Gazetteer; emperor penguin colonies 2023 (Fretwell, 2024, UK Polar " *
        "Data Centre); Adélie penguin colonies from MAPPPD (Humphries et al., 2017), grouped within 100 km. " :
        "Topography from BedMachine Greenland v6 (Morlighem et al., 2017) and ETOPO 2022; ice velocity from the " *
        "MEaSUREs multi-year mosaic (Joughin et al., 2018); median sea-ice edge 1981–2010 from the Sea Ice Index v4 " *
        "(Fetterer et al., 2017); place names from GeoNames; wildlife areas from the Areas Important to Wildlife " *
        "(Greenland Institute of Natural Resources, DCE), one sign per area, grouped within 100 km. "
    return base * "Expedition routes from the expedition accounts and published positions " *
           "(data/routes_$(region).csv); arcs join documented waypoints only. Polar stereographic projection."
end

function main(args=ARGS)
    isempty(args) && error("usage: vintage.jl antarctica|greenland [noroutes] [nofauna] [full] [pdf]")
    region = args[1]
    region in ("antarctica", "greenland") || error("unknown region $region")
    flags = args[2:end]
    for a in flags
        a in ("noroutes", "nofauna", "full", "pdf") || error("unknown option $a")
    end
    full = "full" in flags
    labels = filter(l -> l.tier <= 1, read_labels(joinpath(ROOT, "data", "labels_$(region).csv");
                                                  namecol=region == "greenland" ? "kl" : "name"))
    routes = "noroutes" in flags ? nothing : read_routes(region)
    fauna  = "nofauna" in flags ? nothing : read_fauna(region)
    d = load_prepared(region; full)
    plot = region == "antarctica" ? plot_vintage_antarctica : plot_vintage_greenland
    fig = plot(d; labels, routes, fauna)
    tag = join(filter(!isempty, ["noroutes" in flags ? "noroutes" : "", "nofauna" in flags ? "nofauna" : ""]), "_")
    name = "$(region)_A0_vintage" * (isempty(tag) ? "" : "_" * tag) * (full ? "" : "_coarse")
    save_poster(fig, joinpath(ROOT, "plots", "vintage", name); pdf=(full || "pdf" in flags))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
