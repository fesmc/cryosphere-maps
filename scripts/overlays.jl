# Map overlays: historic expedition routes (data/routes_<region>.csv) and
# wildlife glyphs (data/prepared/fauna_<region>.csv, see prepare_overlays.jl).
# Draw order on a map: routes! (fixed tracks; returns their boxes and site labels),
# then the labels, then fauna! (glyphs move off the label boxes). Both take
# `nothing` for no overlay.

# ---------------------------------------------------------------------------
# Glyphs (SVG paths, y down; Makie fits them to the marker size)
# ---------------------------------------------------------------------------

glyph(path) = BezierPath(path; fit=true, flipy=true)

const G_PENGUIN = glyph("M5,0 C6.6,0 7.4,1.2 7.2,2.8 C8.8,4.5 9.6,8 9.2,11.5 C9,13.5 8,15 7,15.4 " *
    "L8.6,16 L1.4,16 L3,15.4 C2,15 1,13.5 0.8,11.5 C0.4,8 1.2,4.5 2.8,2.8 C2.6,1.2 3.4,0 5,0 Z " *
    "M7.1,1.7 L9.2,2.3 L7.1,2.7 Z")
const G_BELLY   = glyph("M5,0 C7.4,0.5 8,4 7.7,7 C7.4,9.5 6.3,10.4 5,10.4 C3.7,10.4 2.6,9.5 2.3,7 C2,4 2.6,0.5 5,0 Z")
const G_BEAR    = glyph("M2,9 C2,6 5,4 9,4 L19,4 C21,3 23,2.5 25,3 L27,3.5 C28.5,4 29.5,5 29.5,6 L28,6.5 " *
    "C27,7 26,7 25,7.5 L24,9 L24,15 L21.5,15 L21,10.5 L19,11 L18.5,15 L16,15 L16,11 L10,11 L9.5,15 L7,15 " *
    "L7,11 C6,11.5 5.5,12 5,15 L2.5,15 L3,11 C2.2,10.5 2,10 2,9 Z")
const G_MUSKOX  = glyph("M3,8 C3,5 6,3 10,3 C13,2 16,2.5 18,3.5 L21,4 C23,4 24.5,5 25,7 L26,10 C26,11 25,11.5 " *
    "24,11 L22.5,10 L21.5,12 L20.5,15 L19.5,13 L18.5,15.5 L17,13.2 L15.5,15 L14,13.2 L12.5,15 L11,13.2 " *
    "L9.5,15.5 L8,13.2 L6.5,15 L5,12.5 C3.5,11.5 3,10 3,8 Z")
const G_WALRUS  = glyph("M4,6 C4,3.5 6,2 8,2.5 C10,3 11,4 13,5 C18,4 24,6 28,10 L30,9 L30,12 L26,12.5 " *
    "C20,13.5 12,13.5 8,12 C6,11 4,9 4,6 Z M5.5,8 L6,13.5 L6.8,8.2 Z M7.5,8.3 L8,13 L8.6,8.5 Z")
const G_NARWHAL = glyph("M6,6 C6,3.5 10,2.5 16,3 C22,3.5 26,5 28,6 L31,3.5 L30.5,6.5 L31,9.5 L28,7.5 " *
    "C25,9 20,9.5 14,9.3 C9,9 6,8 6,6 Z M6.5,5.6 L0,4.6 L6.5,6.2 Z")
const G_SHIP    = glyph("M1,10 L19,10 L16,14 L4,14 Z M9.5,10 L9.5,0 L10.5,0 L10.5,10 Z " *
    "M10.5,1 L17,8 L10.5,8 Z M9.5,2 L4,8.5 L9.5,8.5 Z")

# species => (glyph, size in pt at scale 1, legend text, aspect width/height of the glyph,
#             opacity of the ink)
const FAUNA = Dict(
    "emperor"   => (G_PENGUIN, 26, "Emperor penguin colony", 9.6/16, 1.0),
    "adelie"    => (G_PENGUIN, 16, "Adélie penguin colonies", 9.6/16, 0.55),
    "polarbear" => (G_BEAR,    26, "Polar bear denning area", 27.5/12, 1.0),
    "muskox"    => (G_MUSKOX,  24, "Musk-ox calving area", 23/13.5, 1.0),
    "walrus"    => (G_WALRUS,  26, "Walrus haul-out", 26/11.5, 1.0),
    "narwhal"   => (G_NARWHAL, 30, "Narwhal summer area", 31/6.5, 1.0),
)
const FAUNA_ORDER = ["emperor", "adelie", "polarbear", "muskox", "walrus", "narwhal"]

"Wildlife points (species, name, lat, lon, weight) of prepare_overlays.jl."
function read_fauna(region)
    d, h = readdlm(joinpath(REPO_PREP_DIR, "fauna_$(region).csv"), ',', Any; header=true, quotes=true)
    c(n) = d[:, findfirst(==(n), vec(h))]
    return [(; species=String(s), name=string(n), lat=Float64(la), lon=Float64(lo), weight=Float64(w))
            for (s, n, la, lo, w) in zip(c("species"), c("name"), c("lat"), c("lon"), c("weight"))]
end

"Glyph (a penguin with a light belly, other animals as silhouettes) at the points `pts`."
function fauna_glyph!(ax, species, pts; size, ink=:black, paper=:white)
    g = FAUNA[species]
    scatter!(ax, pts; marker=g[1], markersize=size, color=(ink, g[5]))
    species in ("emperor", "adelie") &&
        scatter!(ax, pts; marker=G_BELLY, markersize=0.44size, color=paper, marker_offset=Vec2f(0, -0.08size))
end

"""
Wildlife glyphs, the larger species first; each glyph moves to the nearest
free spot if it would overlap one already drawn or an obstacle. Returns the
glyph boxes and legend entries; `names` replaces legend texts by species.
"""
fauna!(ax, ::Nothing, proj; kw...) = (Tuple[], [])

function fauna!(ax, fauna, proj; kmpp, scale=1.0, ink=:black, paper=:white, obstacles=Tuple[],
                names=Dict{String, String}())
    boxes = Tuple[]; legend = []
    for sp in FAUNA_ORDER
        pts = [project(proj, f.lon, f.lat) for f in fauna if f.species == sp]
        isempty(pts) && continue
        g = FAUNA[sp]; size = g[2]*scale
        wh = 1.1size*kmpp .* (g[4] >= 1 ? (1.0, 1/g[4]) : (g[4], 1.0))     # fit=true: longer side = size
        pos = copy(pts)
        place_area!(pos, fill(wh, length(pos)), eachindex(pos), vcat(obstacles, boxes); maxshift=fill(1.5, length(pos)))
        fauna_glyph!(ax, sp, Point2f.(pos); size, ink, paper)
        append!(boxes, [boxat(p, wh, :center) for p in pos])
        push!(legend, (MarkerElement(marker=g[1], color=(ink, g[5]), markersize=min(size, 30)), get(names, sp, g[3])))
    end
    return boxes, legend
end

# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

"""
A route from data/routes_<region>.csv: positions (lat, lon) in travel order
(rows with seq > 0) and named sites (rows with a site name; seq = 0 for sites
off the drawn track). `track` is :camps for transcribed positions (drawn as
straight legs) or :waypoints for documented places only (drawn as arcs, to show
the simplification). `mode` is sledge, ship, drift or air.
"""
struct Route
    name::String
    years::String
    mode::String
    track::Symbol
    lat::Vector{Float64}
    lon::Vector{Float64}
    sites::Vector{Tuple{String, Float64, Float64}}
end

function read_routes(region)
    d, h = readdlm(joinpath(ROOT, "data", "routes_$(region).csv"), ',', Any; header=true, quotes=true)
    c(n) = d[:, findfirst(==(n), vec(h))]
    name, years, mode, track, seq = string.(c("route")), string.(c("years")), string.(c("mode")),
                                    string.(c("track")), Int.(c("seq"))
    lat, lon, site = Float64.(c("lat")), Float64.(c("lon")), strip.(string.(c("site")))
    routes = Route[]
    for n in unique(name)
        J = findall(==(n), name)
        I = sort(filter(i -> seq[i] > 0, J); by=i -> seq[i])
        i = J[1]
        push!(routes, Route(n, years[i], mode[i], Symbol(track[i]), lat[I], lon[I],
                            [(site[k], lat[k], lon[k]) for k in J if !isempty(site[k])]))
    end
    return routes
end

"Projected path of a route: straight legs between camps, or arcs bowing by `arc` of each leg between waypoints."
function route_path(r::Route, proj; arc=0.08, n=24)
    P = [Point2f(project(proj, lo, la)...) for (la, lo) in zip(r.lat, r.lon)]
    length(P) == 1 && return P
    pts = Point2f[]
    for k in 1:length(P)-1
        a, b = P[k], P[k+1]
        m = (a + b)/2 + (r.track === :waypoints ? arc*Point2f(-(b - a)[2], (b - a)[1]) : Point2f(0, 0))
        for t in range(0, 1, length=n+1)[1:end-1]
            push!(pts, (1 - t)^2*a + 2t*(1 - t)*m + t^2*b)      # quadratic Bézier (straight when m is the midpoint)
        end
    end
    push!(pts, P[end])
    return pts
end

"""
Name written along the route, beside the track: the stretch of the route where
the rotated text box hits the fewest obstacles and the track is straightest.
Returns the (axis-aligned) box of the label.
"""
function track_label!(ax, pts, txt; kmpp, fontsize, font, color, obstacles, limits=nothing, gap=6.0)
    w, h = textbox(txt, fontsize, kmpp; font)
    s = vcat(0.0, cumsum([hypot((pts[k+1] - pts[k])...) for k in 1:length(pts)-1]))
    L = s[end]; L < 1.2w && return nothing
    at(sv) = (k = clamp(searchsortedlast(s, sv), 1, length(pts) - 1); t = (sv - s[k])/max(s[k+1] - s[k], 1e-9);
              pts[k] + t*(pts[k+1] - pts[k]))
    best = nothing; bestc = Inf
    for sc in range(0.6w, L - 0.6w, length=41), side in (1, -1)
        a, b = at(sc - w/2), at(sc + w/2)
        θ = atan(b[2] - a[2], b[1] - a[1]); abs(θ) > π/2 && (θ -= sign(θ)*π)
        dir = Point2f(cos(θ), sin(θ)); nrm = Point2f(-dir[2], dir[1])
        # deviation of the track from the chord over the label length (km)
        dev = maximum(abs(sum((at(sv) - a) .* nrm)) for sv in range(sc - w/2, sc + w/2, length=9))
        c = (a + b)/2 + side*(h/2 + gap*kmpp + dev)*nrm
        hw = (abs(w*cos(θ)) + abs(h*sin(θ)))/2; hh = (abs(w*sin(θ)) + abs(h*cos(θ)))/2
        bx = (c[1] - hw, c[1] + hw, c[2] - hh, c[2] + hh)
        cost = sum((overlap_area(bx, o) for o in obstacles); init=0.0) + 2dev*w + 0.02abs(sc - L/2)*h +
               (limits === nothing ? 0.0 : 10*outside_area(bx, limits))
        cost < bestc && ((best, bestc) = ((c, θ, bx), cost))
    end
    c, θ, bx = best
    text!(ax, c[1], c[2]; text=txt, fontsize, font, color, rotation=θ, align=(:center, :center))
    return bx
end

const ROUTE_STYLE = Dict("sledge" => :solid, "ship" => :dash, "drift" => :dot, "air" => :dashdot)

"""
Expedition routes in the colours `colors`, with a dot at the start, a ship at
ship positions and the route names along the tracks. Returns the obstacle
boxes (tracks and names), the site labels (PlaceLabels of type "site", for
draw_labels!) and legend entries.
"""
routes!(ax, ::Nothing, proj; kw...) = (Tuple[], PlaceLabel[], [])

function routes!(ax, routes, proj; kmpp, colors, scale=1.0, font=:italic, fontsize=15, lw=2.0,
                 paper=:white, obstacles=Tuple[], limits=nothing)
    paths = [route_path(r, proj) for r in routes]
    tracks = [polyline_boxes(p; step=15kmpp*scale, hw=3kmpp*scale) for p in paths]
    boxes = Tuple[]; sites = PlaceLabel[]; legend = []
    for (k, (r, p)) in enumerate(zip(routes, paths))
        col = colors[mod1(k, length(colors))]
        ls = ROUTE_STYLE[r.mode]
        if length(p) > 1
            lines!(ax, p; color=(paper, 0.8), linewidth=3lw*scale)
            lines!(ax, p; color=col, linewidth=lw*scale, linestyle=ls)
            scatter!(ax, [p[1]]; markersize=7scale, color=col, strokecolor=paper, strokewidth=0.8scale)
        else
            scatter!(ax, p; marker=G_SHIP, markersize=26scale, color=col)
            push!(boxes, boxat(p[1], (26scale*kmpp, 20scale*kmpp), :center))
        end
        txt = "$(r.name) $(r.years)"
        others = vcat(obstacles, boxes, (tracks[j] for j in eachindex(tracks) if j != k)...)
        bx = length(p) > 1 ? track_label!(ax, p, txt; kmpp, fontsize=fontsize*scale, font, color=col,
                                          obstacles=others, limits) : nothing
        bx === nothing || push!(boxes, bx)
        append!(sites, [PlaceLabel(s, lon, lat, "site", 1, "routes_csv", String[]) for (s, lat, lon) in r.sites])
        push!(legend, (LineElement(color=col, linewidth=lw*scale, linestyle=ls), txt))
    end
    return vcat(boxes, tracks...), sites, legend
end

# route colours on the standard posters: dark and warm, apart from the velocity colours
const ROUTE_COLORS = [colorant"#c0142b", colorant"#e66101", colorant"#5c3310", colorant"#1a1a1a",
                      colorant"#b8860b", colorant"#8c510a", colorant"#d6604d"]

