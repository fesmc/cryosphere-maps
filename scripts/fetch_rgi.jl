# Step 3: fetch Randolph Glacier Inventory 6.0 outlines for the regions that
# fall inside the Greenland poster window beyond the BedMachine domain
# (Arctic Canada North/South, Iceland). Mirror: OGGM, Uni Bremen (no login).
#
# Output: data/external/rgi60/<region>/*.shp
# Usage:  julia --project=. scripts/fetch_rgi.jl

using Downloads

const BASE = "https://cluster.klima.uni-bremen.de/~oggm/rgi/www.glims.org/RGI/rgi60_files"
const REGIONS = ["03_rgi60_ArcticCanadaNorth", "04_rgi60_ArcticCanadaSouth", "06_rgi60_Iceland"]
const OUT = joinpath(@__DIR__, "..", "data", "external", "rgi60")

for reg in REGIONS
    dir = joinpath(OUT, reg)
    if isdir(dir) && any(endswith(".shp"), readdir(dir))
        println("have ", reg); continue
    end
    mkpath(dir)
    zip = joinpath(OUT, reg*".zip")
    println("downloading ", reg, " ...")
    Downloads.download("$(BASE)/$(reg).zip", zip)
    run(`unzip -q -o $zip -d $dir`)
    rm(zip)
end
println("RGI outlines in ", normpath(OUT))
