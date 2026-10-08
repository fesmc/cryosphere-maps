# Step 4: rasterise the RGI 6.0 outlines onto the extended 8 km Greenland grid
# (same grid as data/GRL-8KM-EXT_ETOPO2022.nc) as a glacier-ice area fraction,
# using NSUB x NSUB sub-samples per cell (even-odd point-in-polygon test).
#
# Input:  data/external/rgi60/*/*.shp            (from fetch_rgi.jl)
#         data/GRL-8KM-EXT_ETOPO2022.nc          (from prepare_etopo.jl; grid only)
# Output: data/GRL-8KM-EXT_RGI60.nc              (xc, yc [km], ice_frac [0-1])
# Usage:  julia --project=. scripts/prepare_rgi.jl

include("common.jl")
import Shapefile

const RGIDIR = joinpath(ROOT, "data", "external", "rgi60")
const GRIDF  = joinpath(ROOT, "data", "GRL-8KM-EXT_ETOPO2022.nc")
const DST    = joinpath(ROOT, "data", "GRL-8KM-EXT_RGI60.nc")
const NSUB   = 4

"Even-odd point-in-polygon test over all rings (handles holes)."
function inpoly(px, py, xs, ys, starts)
    inside = false
    for (r, s) in enumerate(starts)
        e = r < length(starts) ? starts[r+1] - 1 : length(xs)
        j = e
        for i in s:e
            if (ys[i] > py) != (ys[j] > py) && px < (xs[j] - xs[i])*(py - ys[i])/(ys[j] - ys[i]) + xs[i]
                inside = !inside
            end
            j = i
        end
    end
    return inside
end

function main()
    x = readaxis(GRIDF, "xc"); y = readaxis(GRIDF, "yc"); dx = x[2] - x[1]
    hits = zeros(Int, length(x), length(y))
    sub = ((1:NSUB) .- 0.5) ./ NSUB .- 0.5          # sub-sample offsets in cell units
    shps = [joinpath(RGIDIR, d, f) for d in readdir(RGIDIR) for f in readdir(joinpath(RGIDIR, d)) if endswith(f, ".shp")]
    for shp in shps
        n = 0
        for g in Shapefile.shapes(Shapefile.Table(shp))
            g === missing && continue
            pr = [project(PROJ_GRL, p.x, p.y) for p in g.points]
            xs = first.(pr); ys = last.(pr)
            (maximum(xs) < x[1] || minimum(xs) > x[end] || maximum(ys) < y[1] || minimum(ys) > y[end]) && continue
            starts = g.parts .+ 1
            i0 = max(1, floor(Int, (minimum(xs) - x[1])/dx)); i1 = min(length(x), ceil(Int, (maximum(xs) - x[1])/dx) + 2)
            j0 = max(1, floor(Int, (minimum(ys) - y[1])/dx)); j1 = min(length(y), ceil(Int, (maximum(ys) - y[1])/dx) + 2)
            for j in j0:j1, i in i0:i1, b in sub, a in sub
                inpoly(x[i] + a*dx, y[j] + b*dx, xs, ys, starts) && (hits[i, j] += 1)
            end
            n += 1
        end
        println(basename(shp), ": ", n, " glaciers in grid")
    end
    frac = min.(hits ./ NSUB^2, 1.0)    # overlapping outlines can double count

    NCDataset(DST, "c") do out
        defVar(out, "xc", x, ("xc",); attrib=["units" => "km"])
        defVar(out, "yc", y, ("yc",); attrib=["units" => "km"])
        defVar(out, "ice_frac", frac, ("xc", "yc"); attrib=["units" => "1", "long_name" => "RGI 6.0 glacier area fraction"])
        out.attrib["source"] = "RGI Consortium (2017): Randolph Glacier Inventory 6.0, doi:10.7265/N5-RGI-60; regions 03, 04, 06"
        out.attrib["method"] = "$(NSUB)x$(NSUB) sub-sampled even-odd point-in-polygon per 8 km cell"
    end
    println("wrote ", DST, "; cells with ice_frac > 0.5: ", count(>(0.5), frac))
end

main()
