# Step 2: regrid the ETOPO 2022 subset onto an extended 8 km polar
# stereographic grid aligned with GRL-8KM (same projection, same nodes, padded
# by NPAD cells on every side). Projected ETOPO points are bin-averaged into
# the 8 km cells.
#
# Input:  data/external/ETOPO2022_60s_greenland_subset.nc   (from fetch_etopo.jl)
# Output: data/GRL-8KM-EXT_ETOPO2022.nc                     (xc, yc [km], z [m])
# Usage:  julia --project=. scripts/prepare_etopo.jl

include("common.jl")

const SRC  = joinpath(ROOT, "data", "external", "ETOPO2022_60s_greenland_subset.nc")
const DST  = joinpath(ROOT, "data", "GRL-8KM-EXT_ETOPO2022.nc")
const GRID = joinpath(ICE_DATA, "Greenland", "GRL-8KM", "GRL-8KM_TOPO-M17.nc")
const NPAD = 50                       # 50 cells = 400 km of padding on each side

function main()
    xg = readaxis(GRID, "xc"); yg = readaxis(GRID, "yc"); dx = xg[2] - xg[1]
    x = collect(xg[1] - NPAD*dx:dx:xg[end] + NPAD*dx)
    y = collect(yg[1] - NPAD*dx:dx:yg[end] + NPAD*dx)

    lon = readaxis(SRC, "lon"); lat = readaxis(SRC, "lat"); z = readvar(SRC, "z")
    s = zeros(length(x), length(y)); n = zeros(Int, length(x), length(y))
    for (j, φ) in enumerate(lat), (i, λ) in enumerate(lon)
        isnan(z[i, j]) && continue
        px, py = project(PROJ_GRL, λ, φ)
        a = round(Int, (px - x[1])/dx) + 1; b = round(Int, (py - y[1])/dx) + 1
        (1 <= a <= length(x) && 1 <= b <= length(y)) || continue
        s[a, b] += z[i, j]; n[a, b] += 1
    end
    zz = [n[i, j] > 0 ? s[i, j]/n[i, j] : NaN for i in eachindex(x), j in eachindex(y)]
    nmiss = count(isnan, zz)
    nmiss > 0 && @warn "$nmiss empty cells (outside the fetched subset?)"

    NCDataset(DST, "c") do out
        defVar(out, "xc", x, ("xc",); attrib=["units" => "km"])
        defVar(out, "yc", y, ("yc",); attrib=["units" => "km"])
        defVar(out, "z", zz, ("xc", "yc"); attrib=["units" => "m", "long_name" => "surface elevation (ETOPO 2022)"])
        out.attrib["grid"] = "GRL-8KM polar stereographic (lat_ts 70N, lon0 45W), padded by $(NPAD) cells"
        out.attrib["method"] = "bin-average of projected ETOPO 2022 60s points (stride 2)"
        out.attrib["source"] = NCDataset(ds -> ds.attrib["reference"], SRC)
    end
    println("wrote ", DST, " (", length(x), " x ", length(y), ")")
end

main()
