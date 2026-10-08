# Data locations and poster grid definitions (shared by fetch, prepare, plot).

# Root for raw downloads and prepared grids (large; not in the repo)
const DATA_DIR = get(ENV, "CRYOMAPS_DATA", "/albedo/work/projects/p_forclima/cryosphere-maps-data")
const RAW_DIR  = joinpath(DATA_DIR, "raw")
const PREP_DIR = joinpath(DATA_DIR, "prepared")
const ROOT     = normpath(joinpath(@__DIR__, ".."))

# Poster grids (polar stereographic, km). `dx` is the raster resolution,
# `dxb` the coarser grid used for drainage-basin divides.
const GRIDS = Dict(
    "greenland"  => (epsg=3413, x=(-1110.0, 1350.0), y=(-3430.0, -590.0), dx=0.5, dxb=4.0),
    "antarctica" => (epsg=3031, x=(-3040.0, 3040.0), y=(-3040.0, 3040.0), dx=1.0, dxb=8.0),
)

prepared_file(region) = joinpath(PREP_DIR, "$(region)_$(round(Int, 1000*GRIDS[region].dx))m.nc")

# Small derived data kept in the repo (gazetteers, sea-ice edges, basin outlines)
# and a copy of the poster grid at every DRAFT_STRIDE-th node, which is enough
# for the 1600 px sharing PNGs and needs no access to $CRYOMAPS_DATA.
const REPO_PREP_DIR = joinpath(ROOT, "data", "prepared")
const DRAFT_STRIDE  = 4
draft_file(region) = joinpath(REPO_PREP_DIR, "$(region)_$(round(Int, 1000*DRAFT_STRIDE*GRIDS[region].dx))m.nc")
