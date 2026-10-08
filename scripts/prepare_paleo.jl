# Paleo step 1: ice sheets at a past time slice from the PaleoMIST 1.0 reconstruction
# (Gowan et al., 2021), as changes relative to the present applied to present-day
# topography: ETOPO 2022 for the Northern Hemisphere, BedMachine Antarctica v4 (the
# repo copy of the poster grid) for Antarctica.
#
# PaleoMIST provides, on 5 km grids in regional projections (North America incl.
# Greenland and Iceland, Eurasia, Antarctica):
#   thickness/<yr>.nc  grounded ice thickness
#   deform/<yr>.nc     change of the bed elevation relative to sea level (deformation
#                      and sea-level change), zero at present
# and, globally on a 0.25° grid, the deformed base topography at every time slice.
#
# With H the ice thickness and zb the bed elevation relative to sea level:
#   zb(t) = zb_pd + deform(t)       (global 0.25° grid outside the regional domains)
#   H(t)  = max(H_pd + H_pm(t) - H_pm(0), 0)  where the base has grounded ice and
#                                             PaleoMIST has ice today,
#           H_pm(t)                           elsewhere,
# so present-day ice keeps the detail of the base data and ice beyond it (and over
# today's ice shelves, which PaleoMIST does not include) comes from PaleoMIST.
# Ice thinner than flotation is floating; PaleoMIST itself has grounded ice only.
#
# Output: data/prepared/paleo_<region>_<t>ka.nc in the repo with z_srf, z_bed [m], mask
# (0 ocean, 1 ice-free land, 2 grounded ice, 3 floating ice) at the time slice and
# pd_mask (the same at present), plus key numbers at the time slice and at present
# (prefix pd_) as global attributes.
#
# Usage: julia --project=. scripts/prepare_paleo.jl [time=20] [nh] [antarctica]
# Needs the raw data of step 0 (scripts/fetch_data.jl); takes a few minutes.

include("prepare.jl")

const PM_DIR = raw("paleomist", "Gowan_ice_reconstruction", "ice_reconstruction")
pm(parts...) = joinpath(PM_DIR, parts...)

# Projections of the regional PaleoMIST grids (GMT, see <region>/projection_info.sh):
# grid coordinates are metres from the lower-left corner (lon, lat) of the GMT region,
# except for Antarctica, which is EPSG:3031 already.
const PM_REGIONS = Dict(
    "North_America" => (srs="+proj=laea +lat_0=60 +lon_0=-94 +ellps=WGS84", corner=(-135.0, 25.0)),
    "Eurasia"       => (srs="+proj=laea +lat_0=90 +lon_0=0 +ellps=WGS84", corner=(-12.0, 47.0)),
    "Antarctica"    => (srs="EPSG:3031", corner=nothing),
)

"Projected coordinates (m) of (lon, lat) in `srs`."
function transform(srs, lon, lat)
    out = read(pipeline(IOBuffer("$lon $lat\n"), `$(gdaltransform_exe()) -s_srs EPSG:4326 -t_srs $srs`), String)
    v = parse.(Float64, split(out))
    return v[1], v[2]
end

"""
Copy of a regional PaleoMIST grid with projected coordinates, which gdalwarp can
read; returns the GDAL source and its projection.
"""
function pm_georef(region, var, file)
    p = PM_REGIONS[region]
    src = pm("ice_reconstruction_files", region, var, file)
    x, y, z = NCDataset(ds -> (ds["x"][:], ds["y"][:], ds["z"][:, :]), src)
    x0, y0 = p.corner === nothing ? (0.0, 0.0) : transform(p.srs, p.corner...)
    mkpath(TMP); out = joinpath(TMP, "pm_$(region)_$(var)_$(file)")
    rm(out; force=true)
    NCDataset(out, "c") do ds
        defVar(ds, "x", x .+ x0, ("x",); attrib=["units" => "m", "standard_name" => "projection_x_coordinate"])
        defVar(ds, "y", y .+ y0, ("y",); attrib=["units" => "m", "standard_name" => "projection_y_coordinate"])
        defVar(ds, "z", Float32.(coalesce.(z, NaN32)), ("x", "y"); attrib=["_FillValue" => NaN32])
    end
    return "NETCDF:$(out):z", p.srs
end

"Regional PaleoMIST field `var` at `yr` (years BP) on grid `g` (NaN outside the region)."
function pm_field(region, var, yr, g; resample="average")
    file = var == "modern_topo" ? "$(region).nc" : "$(yr).nc"
    src, srs = pm_georef(region, var, file)
    return warp(src, g; s_srs=srs, resample, name="pm_$(region)_$(var)_$(yr)")
end

"""
Change of the PaleoMIST base topography (deformation and sea-level change) between
`yr` and the present on the global 0.25° grid, on grid `g`.
"""
function pm_global_deform(yr, g)
    f = pm("global_grid", "reconstruction_0.25_degree.nc")
    lon, lat, d = NCDataset(f) do ds
        t = ds["time"].var[:]                    # years relative to 1950 (unparseable units)
        it, i0 = findfirst(==(-yr), t), findfirst(==(0), t)
        it === nothing && error("no time slice $yr BP in $f")
        b = ds["base_topography"].var
        ds["lon"][:], ds["lat"][:], Float32.(coalesce.(b[:, :, it] .- b[:, :, i0], NaN32))
    end
    keep = lon .< 180                            # drop the repeated 180° column
    mkpath(TMP); out = joinpath(TMP, "pm_global_deform_$(yr).nc"); rm(out; force=true)
    NCDataset(out, "c") do ds
        defVar(ds, "lon", lon[keep], ("lon",); attrib=["units" => "degrees_east", "standard_name" => "longitude"])
        defVar(ds, "lat", lat, ("lat",); attrib=["units" => "degrees_north", "standard_name" => "latitude"])
        defVar(ds, "z", d[keep, :], ("lon", "lat"); attrib=["_FillValue" => NaN32])
    end
    return warp("NETCDF:$(out):z", g; s_srs="EPSG:4326", resample="bilinear", name="pm_global_deform_$(yr)")
end

"Value of the first of the fields that is finite at each cell (NaN if none)."
function firstfinite(fields...)
    out = fill(NaN, size(first(fields)))
    for f in reverse(fields)
        ok = isfinite.(f); out[ok] .= f[ok]
    end
    return out
end

"""
PaleoMIST ice thickness at `yr` and at present, and the bed change at `yr`, on
grid `g` from the regional grids `regions` (the global grid outside them).
"""
function pm_fields(regions, yr, g)
    z0 = zeros(length(axes_km(g)[1]), length(axes_km(g)[2]))
    H, H0 = copy(z0), copy(z0)
    for reg in regions
        H  = max.(H,  replace(pm_field(reg, "thickness", yr, g), NaN => 0.0))
        H0 = max.(H0, replace(pm_field(reg, "thickness", 0, g), NaN => 0.0))
    end
    D = firstfinite([pm_field(reg, "deform", yr, g; resample="bilinear") for reg in regions]...,
                    pm_global_deform(yr, g))
    return H, H0, D
end

"""
Check the georeferencing: RMS difference between the PaleoMIST present-day
topography and the base bed `zb` where both are defined.
"""
function check_registration(regions, g, zb)
    for reg in regions
        mt = pm_field(reg, "modern_topo", 0, g)
        ok = isfinite.(mt) .& isfinite.(zb)
        d = mt[ok] .- zb[ok]
        println("  $reg: PaleoMIST modern topography - base bed: mean $(round(sum(d)/length(d), digits=1)) m, ",
                "rms $(round(sqrt(sum(abs2, d)/length(d)), digits=1)) m over $(count(ok)) cells")
    end
end

const H_MIN = 10.0     # m, thinner ice is left out
const NUMBER_KEYS = ["ice_area_km2", "floating_area_km2", "ice_volume_km3", "sea_level_equivalent_m", "max_thickness_m"]

"""
Ice thickness, bed, surface and mask at the time slice from the base (bed `zb`,
thickness `Hpd`, grounded ice `grounded`) and the PaleoMIST fields.
"""
function apply_paleo(zb, Hpd, grounded, H, H0, D)
    anom = grounded .& (H0 .> 0)
    Ht = ifelse.(anom, max.(Hpd .+ H .- H0, 0.0), H)
    zbt = zb .+ D
    mask = zeros(size(zbt)); zs = similar(zbt)
    for i in eachindex(zbt)
        h, b = Ht[i], zbt[i]
        if h < H_MIN
            Ht[i] = 0.0
            mask[i] = b < 0 ? 0 : 1
            zs[i] = max(b, 0.0)
        elseif h*RHO_ICE/RHO_SEA < -b
            mask[i] = 3; zs[i] = h*(1 - RHO_ICE/RHO_SEA)
        else
            mask[i] = 2; zs[i] = b + h
        end
    end
    return Ht, zbt, zs, mask
end

"Key numbers of the ice in each (prefix, domain) of `domains`, prefixed with `pre`."
region_numbers(g, proj, domains, mask, H, zb; pre="") =
    reduce(vcat, [ice_numbers(g, proj, mask, H, zb; domain=dom, prefix=pre*p) for (p, dom) in domains])

"Change in sea-level equivalent from the present (prefix pd_) to the time slice, for each region prefix."
sle_change(numbers, prefixes) = (d = Dict(numbers);
    [p*"sle_change_m" => d[p*"sea_level_equivalent_m"] - d["pd_"*p*"sea_level_equivalent_m"] for p in prefixes])

"Is point p inside polygon `poly` (vectors of (x, y))? Even-odd rule."
function inpoly(p, poly)
    c = false; n = length(poly)
    for i in 1:n
        a, b = poly[i], poly[mod1(i + 1, n)]
        ((a[2] > p[2]) != (b[2] > p[2])) &&
            p[1] < (b[1] - a[1])*(p[2] - a[2])/(b[2] - a[2]) + a[1] && (c = !c)
    end
    return c
end

# Approximate outline (lon, lat) separating Greenland from Iceland and the ice of
# Arctic Canada (through Davis Strait, Baffin Bay and Nares Strait)
const GREENLAND_POLY = [(-44, 58.5), (-54, 62), (-58.5, 66.5), (-64, 72), (-73, 76.5), (-73.8, 78.3), (-69.5, 79.7),
                        (-66.5, 80.8), (-62.5, 81.6), (-59.5, 82.3), (-50, 84), (-20, 85), (-3, 80), (-8, 72),
                        (-29, 66.5), (-38, 60)]

"""
Region of every cell of the NH grid: Greenland (GREENLAND_POLY), Iceland (a box
around it), Eurasia (12°W–150°E) and North America (the rest).
"""
function nh_regions(g)
    x, y = axes_km(g)
    poly = [project(PROJ_GRL, lon, lat) for (lon, lat) in GREENLAND_POLY]
    reg = fill("", length(x), length(y))
    for j in eachindex(y), i in eachindex(x)
        lon, lat = unproject(PROJ_GRL, x[i], y[j])
        lon = mod(lon + 180, 360) - 180
        reg[i, j] = inpoly((x[i], y[j]), poly) ? "greenland" :
                    (-26 <= lon <= -12 && 62.5 <= lat <= 67.5) ? "iceland" :
                    (-12 <= lon < 150) ? "eurasia" : "namerica"
    end
    return reg
end

"Northern Hemisphere: PaleoMIST North America and Eurasia on ETOPO 2022."
function prepare_nh(t)
    g = PALEO_GRIDS["nh"]; yr = round(Int, 1000t)
    et(v) = warp(bm_var(raw("etopo", "ETOPO2022_60s_$(v)_nh.nc"), "z"), g; s_srs="EPSG:4326", name="nh_etopo_$v")
    zs_pd, zb_pd = et("surface"), et("bed")
    regions = ["North_America", "Eurasia"]
    check_registration(regions, g, zb_pd)
    H, H0, D = pm_fields(regions, yr, g)

    # present: ETOPO ice (surface above bed: Greenland) and the PaleoMIST ice caps
    Hpd = max.(zs_pd .- zb_pd, 0.0)
    grounded = Hpd .>= H_MIN
    pdmask = ifelse.(grounded .| (H0 .>= H_MIN), 2.0, ifelse.(zs_pd .< 0, 0.0, 1.0))
    Ht, zbt, zst, mask = apply_paleo(zb_pd, Hpd, grounded, H, H0, D)

    reg = nh_regions(g)
    domains = [("", trues(size(reg))); [(r*"_", reg .== r) for r in ("namerica", "greenland", "iceland", "eurasia")]]
    numbers = vcat(region_numbers(g, PROJ_GRL, domains, mask, Ht, zbt),
                   region_numbers(g, PROJ_GRL, domains, pdmask, Hpd, zb_pd; pre="pd_"))
    append!(numbers, sle_change(numbers, first.(domains)))
    write_paleo("nh", t, g, zst, zbt, mask, pdmask, vcat(numbers,
        ["sources" => "PaleoMIST 1.0 (North_America, Eurasia; global 0.25 degree grid for the bed change elsewhere); " *
                      "ETOPO 2022 60s surface and bedrock elevation",
         "regions" => "namerica, greenland, iceland, eurasia: approximate, see nh_regions in scripts/prepare_paleo.jl"]))
end

"Antarctica: PaleoMIST Antarctica on BedMachine Antarctica v4 (repo copy of the poster grid)."
function prepare_antarctica_paleo(t)
    g = PALEO_GRIDS["antarctica"]; yr = round(Int, 1000t)
    f = draft_file("antarctica")
    x, y, zs_pd, zb_pd, pdmask, pd = NCDataset(f) do ds
        rd(v) = Float64.(coalesce.(ds[v][:, :], NaN))
        ds["x"][:], ds["y"][:], rd("z_srf"), rd("z_bed"), rd("mask"), Dict(ds.attrib)
    end
    (x[1], x[end], y[1], y[end], x[2] - x[1]) == (g.x..., g.y..., g.dx) || error("PALEO_GRIDS[\"antarctica\"] does not match $f")
    Hpd = ifelse.(pdmask .== 2, max.(zs_pd .- zb_pd, 0.0), 0.0)     # grounded ice; floating ice is not needed
    check_registration(["Antarctica"], g, zb_pd)
    H, H0, D = pm_fields(["Antarctica"], yr, g)
    Ht, zbt, zst, mask = apply_paleo(zb_pd, Hpd, pdmask .== 2, H, H0, D)

    # East, West and the Peninsula as for the present-day poster, with the ice at the
    # time slice (and today's) given to the nearest region
    names = ["East", "West", "Peninsula"]
    reg = imbie_regions(g, "ANT", "Regions", names, (mask .>= 2) .| (pdmask .>= 2); detached=true, name="ant_regions")
    domains = [("", trues(size(reg))); [(lowercase(n)*"_", reg .== k) for (k, n) in enumerate(names)]]
    # present: the key numbers of the present-day poster (full-resolution grid, BedMachine thickness)
    numbers = vcat(region_numbers(g, PROJ_ANT, domains, mask, Ht, zbt),
                   ["pd_"*k => v for (k, v) in pd if any(endswith(k, s) for s in NUMBER_KEYS)])
    append!(numbers, sle_change(numbers, first.(domains)))
    write_paleo("antarctica", t, g, zst, zbt, mask, pdmask, vcat(numbers,
        ["sources" => "PaleoMIST 1.0 (Antarctica); BedMachine Antarctica v4 ($(basename(f)))",
         "regions" => "east, west, peninsula: IMBIE 2 regions, ice outside them given to the nearest region"]))
end

function write_paleo(region, t, g, zs, zb, mask, pdmask, attrs)
    x, y = axes_km(g)
    write_grids(paleo_file(region, t), x, y,
        ["z_srf" => (zs, "m"), "z_bed" => (zb, "m"), "mask" => (mask, "1"), "pd_mask" => (pdmask, "1")],
        vcat(attrs, ["time_ka" => Float64(t), "grid" => "EPSG:$(g.epsg), $(g.dx) km",
                     "method" => "PaleoMIST 1.0 changes relative to the present applied to present-day topography " *
                                 "(see scripts/prepare_paleo.jl)"]))
end

function main_paleo(args)
    i = findfirst(startswith("time="), args)
    t = i === nothing ? 20.0 : parse(Float64, split(args[i], "=")[2])
    regions = filter(!startswith("time="), args)
    for reg in (isempty(regions) ? ["nh", "antarctica"] : regions)
        println("== ", reg, " at ", katag(t))
        reg == "nh"         ? prepare_nh(t) :
        reg == "antarctica" ? prepare_antarctica_paleo(t) : error("unknown region $reg")
    end
    rm(TMP; recursive=true, force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_paleo(ARGS)
end
