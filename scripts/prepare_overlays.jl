# Wildlife points for the map overlays (scripts/overlays.jl), from the raw files of
# `scripts/fetch_data.jl fauna`, to data/prepared/fauna_{antarctica,greenland}.csv:
#   - emperor penguin colonies 2023 (BAS), one point each;
#   - Adélie penguin colonies (MAPPPD), grouped within CLUSTER_KM["adelie"];
#   - Greenland Areas Important to Wildlife (GINR/DCE), one point per area at its
#     centroid, grouped within CLUSTER_KM.
# Usage: julia --project=. scripts/prepare_overlays.jl

import JSON
using DelimitedFiles, Printf
include("paths.jl")
include("projection.jl")

const FAUNA_RAW = joinpath(RAW_DIR, "fauna")
const CLUSTER_KM = Dict("adelie" => 100.0, "walrus" => 100.0, "narwhal" => 100.0, "polarbear" => 100.0,
                        "muskox" => 100.0)
const AIW_LAYERS = [("muskox", "muskox_calving"), ("walrus", "walrus_haulout"),
                    ("polarbear", "polarbear_denning"), ("narwhal", "narwhal_summer")]

# a point: species, name, projected (x, y) in km, weight (count or area), source
const Pt = @NamedTuple{species::String, name::String, xy::Tuple{Float64, Float64}, w::Float64, source::String}

"Emperor penguin colonies (name, lon, lat) from the BAS KMZ."
function emperor_colonies()
    kml = read(`unzip -p $(joinpath(FAUNA_RAW, "emperor_colony_locations2023.kmz")) doc.kml`, String)
    pts = Pt[]
    for m in eachmatch(r"<Placemark.*?<name>(.*?)</name>.*?<coordinates>\s*([-\d.]+),([-\d.]+)"s, kml)
        lon, lat = parse(Float64, m[2]), parse(Float64, m[3])
        push!(pts, (; species="emperor", name=String(m[1]), xy=project(PROJ_ANT, lon, lat), w=1.0,
                      source="Fretwell (2024), BAS"))
    end
    return pts
end

"Adélie penguin sites from MAPPPD, weighted by their latest count."
function adelie_sites()
    d, h = readdlm(joinpath(FAUNA_RAW, "mapppd_AllCounts.csv"), ',', Any; header=true, quotes=true)
    c(n) = findfirst(==(n), vec(h))
    latest = Dict{String, Tuple{Int, Pt}}()
    for r in eachrow(d)
        r[c("common_name")] == "adelie penguin" || continue
        id = string(r[c("site_id")]); yr = Int(r[c("year")])
        n = r[c("penguin_count")]; n = n isa Number ? Float64(n) : 1.0
        p = (; species="adelie", name=string(r[c("site_name")]),
               xy=project(PROJ_ANT, Float64(r[c("longitude_epsg_4326")]), Float64(r[c("latitude_epsg_4326")])),
               w=n, source="MAPPPD")
        (!haskey(latest, id) || yr > latest[id][1]) && (latest[id] = (yr, p))
    end
    return [p for (_, p) in values(latest)]
end

"Area-weighted centroid (km) and area (km²) of the polygons of a GeoJSON geometry."
function centroid(g, proj)
    polys = g["type"] == "Polygon" ? [g["coordinates"]] : g["coordinates"]
    A = 0.0; cx = 0.0; cy = 0.0
    for poly in polys
        ring = [project(proj, Float64(c[1]), Float64(c[2])) for c in poly[1]]   # outer ring only
        for k in 1:length(ring)-1
            (x0, y0), (x1, y1) = ring[k], ring[k+1]
            a = x0*y1 - x1*y0
            A += a/2; cx += (x0 + x1)*a/6; cy += (y0 + y1)*a/6
        end
    end
    return (cx/A, cy/A), abs(A)
end

"One point per Area Important to Wildlife, weighted by its area."
function aiw_areas(species, layer)
    pts = Pt[]
    for f in JSON.parsefile(joinpath(FAUNA_RAW, "aiw_$(layer).geojson"))["features"]
        f["geometry"] === nothing && continue
        xy, a = centroid(f["geometry"], PROJ_GRL)
        push!(pts, (; species, name=something(f["properties"]["GeographicalName"], ""), xy, w=a,
                      source="Areas Important to Wildlife (GINR/DCE)"))
    end
    return pts
end

"""
Group points closer than `dkm`: the heaviest point (whose name the group keeps)
takes the weight of the others within `dkm` and moves to their weighted mean
position.
"""
function cluster(pts, dkm)
    left = sort(pts; by=p -> -p.w); out = Pt[]
    while !isempty(left)
        p = popfirst!(left)
        near = [q for q in left if hypot((q.xy .- p.xy)...) < dkm]
        filter!(q -> !(q in near), left)
        grp = vcat(p, near); W = sum(q.w for q in grp)
        xy = (sum(q.w*q.xy[1] for q in grp)/W, sum(q.w*q.xy[2] for q in grp)/W)
        push!(out, (; p.species, p.name, xy, w=W, p.source))
    end
    return out
end

function write_fauna(region, pts, proj)
    out = joinpath(REPO_PREP_DIR, "fauna_$(region).csv")
    open(out, "w") do io
        println(io, "species,name,lat,lon,weight,source")
        for p in sort(pts; by=p -> (p.species, -p.w))
            lon, lat = unproject(proj, p.xy...)
            @printf(io, "%s,\"%s\",%.4f,%.4f,%.0f,%s\n", p.species, replace(p.name, "\"" => "'"), lat, lon, p.w, p.source)
        end
    end
    println("wrote ", out, " (", length(pts), " points)")
end

function main()
    ant = vcat(emperor_colonies(), cluster(adelie_sites(), CLUSTER_KM["adelie"]))
    write_fauna("antarctica", ant, PROJ_ANT)
    grl = vcat([cluster(aiw_areas(sp, layer), CLUSTER_KM[sp]) for (sp, layer) in AIW_LAYERS]...)
    write_fauna("greenland", grl, PROJ_GRL)
end

main()
