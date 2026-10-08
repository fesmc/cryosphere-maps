# Step 1: regrid the raw data (step 0) onto the poster grids defined in
# scripts/paths.jl, using GDAL (gdalwarp / ogr2ogr / gdal_rasterize from GDAL_jll).
#
# Output: $CRYOMAPS_DATA/prepared/<region>_<res>m.nc with, on the poster grid,
#   z_srf, z_bed [m], u [m/yr], mask (0 ocean, 1 ice-free land, 2 grounded ice, 3 floating ice)
# and, on the coarser basin grid (xb, yb), basin (integer id, 0 = none; names in
# the `names` attribute). Key numbers of the ice sheet (area, volume, sea-level
# equivalent) are stored as global attributes.
#
# Small derived files go to data/prepared/ in the repo:
#   <region>_<4 x res>m.nc         the same grid at every 4th node (for quick plots)
#   gazetteer_<region>.csv         ice-feature names from GeoNames / SCAR CGA
#   seaice_<region>_<MM>.geojson   median sea-ice edge 1981-2010 for month MM
#
# Usage: julia --project=. scripts/prepare.jl [greenland] [antarctica]
# Takes ~2 min and <3 GB on albedo (jobs/prepare.sh).

using NCDatasets, GDAL_jll, PROJ_jll, DelimitedFiles
include("paths.jl")
include("projection.jl")

const TMP = joinpath(PREP_DIR, "tmp")

# GDAL needs its data files and PROJ's database
ENV["GDAL_DATA"] = joinpath(GDAL_jll.artifact_dir, "share", "gdal")
ENV["PROJ_DATA"] = joinpath(PROJ_jll.artifact_dir, "share", "proj")
ENV["PROJ_LIB"]  = ENV["PROJ_DATA"]

raw(parts...) = joinpath(RAW_DIR, parts...)

"Grid extent in metres as (xmin, ymin, xmax, ymax) of the cell edges."
extent(g; dx=g.dx) = 1e3 .* (g.x[1] - dx/2, g.y[1] - dx/2, g.x[2] + dx/2, g.y[2] + dx/2)

"""
Warp `src` (any GDAL dataset string) onto grid `g` and return the field as a
Float64 matrix (x, y) with south-to-north y, NaN outside the source.
"""
function warp(src, g; s_srs, resample="average", dx=g.dx, name=basename(string(src)))
    mkpath(TMP)
    out = joinpath(TMP, "$(name).nc")
    te = extent(g; dx)
    run(`$(gdalwarp_exe()) -q -overwrite -s_srs $s_srs -t_srs EPSG:$(g.epsg) -te $(te[1]) $(te[2]) $(te[3]) $(te[4])
         -tr $(1e3dx) $(1e3dx) -r $resample -ot Float32 -dstnodata nan -wm 2048 -multi -of netCDF $src $out`)
    return read_gdal_nc(out)
end

"Read Band1 of a GDAL-written netCDF file as (x, y) with ascending y."
function read_gdal_nc(f)
    NCDataset(f) do ds
        z = Float64.(coalesce.(ds["Band1"][:, :], NaN))
        y = ds["y"][:]
        return y[1] > y[end] ? reverse(z, dims=2) : z
    end
end

"""
Rasterise polygons of vector file(s) `src` onto grid `g`, burning value `burn`;
`where` is an optional attribute filter (OGR SQL).
"""
function rasterize(src, g; dx=g.dx, burn=1, name="rasterized", where=nothing)
    mkpath(TMP)
    reproj = joinpath(TMP, "$(name).gpkg")
    out = joinpath(TMP, "$(name).nc")
    srcs = src isa AbstractVector ? src : [src]
    rm(reproj; force=true)
    filt = where === nothing ? `` : `-where $where`
    for s in srcs
        run(`$(ogr2ogr_exe()) -q -t_srs EPSG:$(g.epsg) $filt -append -nln layer -nlt PROMOTE_TO_MULTI $reproj $s`)
    end
    te = extent(g; dx)
    run(`$(gdal_rasterize_exe()) -q -burn $burn -l layer -te $(te[1]) $(te[2]) $(te[3]) $(te[4]) -tr $(1e3dx) $(1e3dx)
         -ot Int16 -init 0 -of netCDF $reproj $out`)
    return replace(read_gdal_nc(out), NaN => 0.0)     # the init value reads back as missing
end

"""
Rasterise named drainage basins: each distinct value of string attribute `attr`
(except `exclude`) gets id 1, 2, ... in sorted order. Returns (ids, names).
"""
function rasterize_basins(shp, g, attr; exclude=String[], name="basins")
    layer = splitext(basename(shp))[1]
    info = read(`$(ogrinfo_exe()) -ro -q -sql "SELECT DISTINCT $attr FROM \"$layer\"" $shp`, String)
    names = sort([strip(m[1]) for m in eachmatch(Regex("$attr \\(String\\) = (.+)"), info)])
    names = filter(n -> n ∉ exclude && n != "(null)", names)
    ids = nothing
    for (k, n) in enumerate(names)
        mkpath(TMP); one = joinpath(TMP, "$(name)_$(k)_src.gpkg")
        run(`$(ogr2ogr_exe()) -q -overwrite -where "$attr = '$n'" $one $shp`)
        r = rasterize(one, g; dx=g.dxb, burn=k, name="$(name)_$k")
        ids = ids === nothing ? r : max.(ids, r)
    end
    return Int.(ids), names
end

axes_km(g; dx=g.dx) = (collect(g.x[1]:dx:g.x[2]), collect(g.y[1]:dx:g.y[2]))

"""
Write the poster grid of `region` to \$CRYOMAPS_DATA/prepared and a copy at every
DRAFT_STRIDE-th node to data/prepared in the repo.
"""
function write_prepared(region, g, fields, basin, attrs)
    x, y = axes_km(g); xb, yb = axes_km(g; dx=g.dxb)
    write_grids(prepared_file(region), x, y, fields,
                vcat(attrs, ["grid" => "EPSG:$(g.epsg), $(g.dx) km (basins $(g.dxb) km)"]); basins=(xb, yb, basin...))
    s = DRAFT_STRIDE
    write_grids(draft_file(region), x[1:s:end], y[1:s:end],
                [k => (v[1:s:end, 1:s:end], u) for (k, (v, u)) in fields],
                vcat(attrs, ["grid" => "EPSG:$(g.epsg), $(s*g.dx) km (basins $(g.dxb) km)",
                             "draft" => "every $(s)th node of $(basename(prepared_file(region)))"]);
                basins=(xb, yb, basin...))
end

"""
Write `fields` (name => (matrix, units); masks, named `*mask`, as Int8) on the
x, y grid (km) and the global attributes `attrs` to `out`; `basins` = (xb, yb,
ids, names) adds the basin ids on their coarser grid.
"""
function write_grids(out, x, y, fields, attrs; basins=nothing)
    mkpath(dirname(out)); rm(out; force=true)
    NCDataset(out, "c") do ds
        defVar(ds, "x", x, ("x",); attrib=["units" => "km"])
        defVar(ds, "y", y, ("y",); attrib=["units" => "km"])
        for (k, (v, units)) in fields
            T = endswith(k, "mask") ? Int8 : Float32
            defVar(ds, k, T.(replace(v, NaN => (T == Int8 ? -1 : NaN32))), ("x", "y"); deflatelevel=4,
                   attrib=["units" => units])
        end
        if basins !== nothing
            xb, yb, basin, basin_names = basins
            defVar(ds, "xb", xb, ("xb",); attrib=["units" => "km"])
            defVar(ds, "yb", yb, ("yb",); attrib=["units" => "km"])
            defVar(ds, "basin", Int16.(basin), ("xb", "yb"); deflatelevel=4,
                   attrib=["units" => "1", "names" => join(basin_names, ",")])
        end
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    println("wrote ", out)
end

# ---------------------------------------------------------------------------
# Key numbers
# ---------------------------------------------------------------------------

"""
Fill the zero cells of an integer label grid with the label of the nearest
labelled cell (multi-source breadth-first search over 8 neighbours), moving
only through cells where `through` is true.
"""
function nearest_fill(lab; through=trues(size(lab)))
    out = copy(lab); nx, ny = size(out)
    q = Tuple{Int32, Int32}[(i, j) for j in 1:ny for i in 1:nx if out[i, j] != 0]
    head = 1
    while head <= length(q)
        i, j = q[head]; head += 1
        for di in -1:1, dj in -1:1
            a, b = i + di, j + dj
            (1 <= a <= nx && 1 <= b <= ny && out[a, b] == 0 && through[a, b]) || continue
            out[a, b] = out[i, j]; push!(q, (a, b))
        end
    end
    return out
end

"""
Region id (1, 2, ... for the values `names` of attribute `attr`; with `others`,
all remaining polygons get id length(names)+1) of every cell of grid `g`,
from the IMBIE 2 polygons of `reg` ("GRE", "ANT"). The polygons cover the
grounded ice, so ice outside them (ice shelves, the ice margin) gets the region
of the nearest polygon it is connected to through ice (`ice`). With `detached`,
ice not connected to any polygon (e.g. islands) gets the nearest region too;
otherwise it is left at 0.
"""
function imbie_regions(g, reg, attr, names, ice; others=false, detached=false, name="regions")
    shp = imbie_shapefile(reg)
    x, y = axes_km(g)
    lab = zeros(Int32, length(x), length(y))
    sel = others ? [names; "others"] : names
    for (k, n) in enumerate(sel)
        where = n == "others" ? join(["$attr <> '$v'" for v in names], " AND ") : "$attr = '$n'"
        r = rasterize(shp, g; where, name="$(name)_$k") .> 0
        lab[r .& (lab .== 0)] .= k
    end
    lab = nearest_fill(lab; through=ice)
    return detached ? nearest_fill(lab) : lab
end

const RHO_ICE    = 917.0       # kg/m³
const RHO_SEA    = 1027.0      # kg/m³, for flotation
const RHO_FRESH  = 1000.0      # kg/m³, melt water volume added to the ocean
const OCEAN_AREA = 3.625e14    # m²

"""
Key numbers of the ice (mask 2 or 3) inside `domain`: area, floating area,
volume, maximum thickness and the sea-level equivalent of the grounded ice
above flotation. Cell areas are corrected for the scale distortion of `proj`.
"""
function ice_numbers(g, proj, mask, H, zb; domain=trues(size(mask)), prefix="")
    x, y = axes_km(g); a0 = (1e3*g.dx)^2
    area = afl = vol = vaf = hmax = 0.0
    for j in eachindex(y), i in eachindex(x)
        (domain[i, j] && (mask[i, j] == 2 || mask[i, j] == 3)) || continue
        a = a0 / scale_factor(proj, x[i], y[j])^2
        h = isnan(H[i, j]) ? 0.0 : H[i, j]
        area += a; vol += a*h; hmax = max(hmax, h)
        if mask[i, j] == 3
            afl += a
        elseif !isnan(zb[i, j])
            vaf += a*max(0.0, h - max(0.0, -zb[i, j])*RHO_SEA/RHO_ICE)
        end
    end
    sle = vaf*RHO_ICE/(RHO_FRESH*OCEAN_AREA)
    imax = argmax(ifelse.(domain .& ((mask .== 2) .| (mask .== 3)) .& isfinite.(H), H, -Inf))
    lon, lat = unproject(proj, x[imax[1]], y[imax[2]])
    println("  $(prefix)max thickness at $(round(lat, digits=2))°, $(round(lon, digits=2))°")
    println("  $(prefix)area $(round(area/1e12, digits=3)) Mkm², volume $(round(vol/1e15, digits=3)) Mkm³, ",
            "SLE $(round(sle, digits=2)) m, max thickness $(round(Int, hmax)) m")
    out = [prefix*"ice_area_km2" => area/1e6, prefix*"floating_area_km2" => afl/1e6, prefix*"ice_volume_km3" => vol/1e9,
           prefix*"sea_level_equivalent_m" => sle, prefix*"max_thickness_m" => hmax]
    isempty(prefix) || return out
    return [out;
            "numbers_note" => "on the poster grid; SLE of grounded ice above flotation with " *
                              "rho_ice = $RHO_ICE, rho_sea = $RHO_SEA (flotation), rho_fresh = $RHO_FRESH kg/m3, " *
                              "ocean area $(OCEAN_AREA/1e6) km2"]
end

speed(vx, vy) = hypot.(vx, vy)

bm_var(file, v) = "NETCDF:$(file):$(v)"

"""
Greenland: BedMachine v6 inside its domain; outside it ETOPO 2022 for topography
and RGI 6.0 for glacier ice. RGI also marks ice caps that BedMachine classifies
as ice-free land (e.g. Ellesmere). Velocity from MEaSUREs (NSIDC-0670).
"""
function prepare_greenland()
    g  = GRIDS["greenland"]; s = "EPSG:3413"
    bm = raw("bedmachine", "BedMachineGreenland-v6.nc")
    zs = warp(bm_var(bm, "surface"), g; s_srs=s, name="grl_surface")
    zb = warp(bm_var(bm, "bed"), g; s_srs=s, name="grl_bed")
    H  = warp(bm_var(bm, "thickness"), g; s_srs=s, name="grl_thickness")
    m  = warp(bm_var(bm, "mask"), g; s_srs=s, resample="mode", name="grl_mask")
    vx = warp(raw("measures", "greenland_vel_mosaic250_vx_v1.tif"), g; s_srs=s, name="grl_vx")
    vy = warp(raw("measures", "greenland_vel_mosaic250_vy_v1.tif"), g; s_srs=s, name="grl_vy")
    et = warp(bm_var(raw("etopo", "ETOPO2022_30s_greenland_subset.nc"), "z"), g;
              s_srs="EPSG:4326", resample="bilinear", name="grl_etopo")
    rgi = rasterize(rgi_shapefiles(), g; name="grl_rgi")

    # key numbers for the ice sheet proper: IMBIE regions without the peripheral ice caps,
    # with ice outside the polygons (ice shelves, margins) given to the nearest connected region
    reg = imbie_regions(g, "GRE", "SUBREGION1", ["ICE_CAP"], (m .== 2) .| (m .== 3); others=true, name="grl_regions")
    numbers = ice_numbers(g, PROJ_GRL, m, H, zb; domain=(reg .== 2))

    out = isnan.(m)                                   # outside the BedMachine domain
    zb[out] .= et[out]; zs[out] .= max.(et[out], 0.0)
    m[out] .= ifelse.(et[out] .< 0, 0.0, 1.0)
    m[(m .== 1) .& (rgi .> 0)] .= 2                   # ice caps outside the ice sheet

    basin = rasterize_basins(imbie_shapefile("GRE"), g, "SUBREGION1"; exclude=["ICE_CAP"], name="grl_basin")
    write_prepared("greenland", g,
        ["z_srf" => (zs, "m"), "z_bed" => (zb, "m"), "u" => (speed(vx, vy), "m/yr"), "mask" => (m, "1")], basin,
        vcat(["sources" => "BedMachine Greenland v6; MEaSUREs NSIDC-0670 v1; ETOPO 2022 30s; RGI 6.0; IMBIE2 basins",
              "numbers_domain" => "grounded and floating ice of the IMBIE 2 regions without ICE_CAP, ice outside " *
                                  "the regions given to the nearest region connected through ice"], numbers))
    seaice_geojson("greenland", g, "N", ("03", "09"))
    gazetteer_greenland(g)
end

"Antarctica: BedMachine v4, MEaSUREs velocity v2 (NSIDC-0484), IMBIE2 basins."
function prepare_antarctica()
    g  = GRIDS["antarctica"]; s = "EPSG:3031"
    bm = raw("bedmachine", "NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc")
    zs = warp(bm_var(bm, "surface"), g; s_srs=s, name="ant_surface")
    zb = warp(bm_var(bm, "bed"), g; s_srs=s, name="ant_bed")
    H  = warp(bm_var(bm, "thickness"), g; s_srs=s, name="ant_thickness")
    m  = warp(bm_var(bm, "mask"), g; s_srs=s, resample="mode", name="ant_mask")
    m[m .== 4] .= 2                                   # Lake Vostok -> grounded ice
    m[isnan.(m)] .= 0
    numbers = ice_numbers(g, PROJ_ANT, m, H, zb)
    # East and West Antarctica and the Peninsula (IMBIE 2 regions; ice shelves and islands
    # are given to the nearest region, through connected ice where possible)
    regions = ["East", "West", "Peninsula"]
    reg = imbie_regions(g, "ANT", "Regions", regions, (m .== 2) .| (m .== 3); detached=true, name="ant_regions")
    for (k, r) in enumerate(regions)
        append!(numbers, ice_numbers(g, PROJ_ANT, m, H, zb; domain=(reg .== k), prefix=lowercase(r)*"_"))
    end
    vel = raw("measures", "antarctica_ice_velocity_450m_v2.nc")
    vx = warp(bm_var(vel, "VX"), g; s_srs=s, name="ant_vx")
    vy = warp(bm_var(vel, "VY"), g; s_srs=s, name="ant_vy")

    basin = rasterize_basins(imbie_shapefile("ANT"), g, "Subregion"; name="ant_basin")
    write_prepared("antarctica", g,
        ["z_srf" => (zs, "m"), "z_bed" => (zb, "m"), "u" => (speed(vx, vy), "m/yr"), "mask" => (m, "1")], basin,
        vcat(["sources" => "BedMachine Antarctica v4; MEaSUREs NSIDC-0484 v2; IMBIE2 basins",
              "numbers_domain" => "all grounded and floating ice in BedMachine; east_, west_, peninsula_: " *
                                  "IMBIE 2 regions, ice outside them (shelves, islands) given to the nearest region"],
             numbers))
    seaice_geojson("antarctica", g, "S", ("02", "09"))
    gazetteer_antarctica(g)
end

# ---------------------------------------------------------------------------
# Sea-ice edges and gazetteers (small; written to data/prepared in the repo)
# ---------------------------------------------------------------------------

"Median sea-ice edge (1981-2010) for the given months, reprojected and clipped to grid `g`."
function seaice_geojson(region, g, hemi, months)
    for mm in months
        dir = raw("seaice", "median_extent_$(hemi)_$(mm)_1981-2010_polyline_v4.0")
        shp = only(filter(endswith(".shp"), readdir(dir; join=true)))
        out = joinpath(REPO_PREP_DIR, "seaice_$(region)_$(mm).geojson")
        mkpath(dirname(out)); rm(out; force=true)
        te = extent(g)
        run(`$(ogr2ogr_exe()) -q -f GeoJSON -t_srs EPSG:$(g.epsg) -clipdst $(te[1]) $(te[2]) $(te[3]) $(te[4])
             -simplify 2000 -lco COORDINATE_PRECISION=0 $out $shp`)
        println("wrote ", out)
    end
end

"Is (lon, lat) inside grid `g`?"
inside(g, proj, lon, lat) = ((x, y) = project(proj, lon, lat); g.x[1] <= x <= g.x[2] && g.y[1] <= y <= g.y[2])

"""
Write a gazetteer CSV with columns name, lat, lon, type, source, tier (3 =
gazetteer), alt (other names, separated by |).
"""
function write_gazetteer(region, rows)
    out = joinpath(REPO_PREP_DIR, "gazetteer_$(region).csv")
    mkpath(dirname(out))
    sort!(rows; by=r -> r[1])
    q(s) = occursin(r"[,\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
    open(out, "w") do io
        println(io, "name,lat,lon,type,source,tier,alt")
        for (name, lat, lon, type, src, alt) in rows
            println(io, join([q(name), round(lat, digits=4), round(lon, digits=4), type, q(src), 3, q(alt)], ","))
        end
    end
    println("wrote ", out, " (", length(rows), " names)")
end

"Greenland glaciers from GeoNames (feature code GLCR)."
function gazetteer_greenland(g)
    d = readdlm(raw("geonames", "GL", "GL.txt"), '\t', Any; quotes=false)
    rows = []
    for r in eachrow(d)
        r[8] == "GLCR" || continue
        lat, lon = Float64(r[5]), Float64(r[6])
        inside(g, PROJ_GRL, lon, lat) || continue
        alt = filter(a -> a != r[2], unique(split(string(r[4]), ",")))
        push!(rows, (string(r[2]), lat, lon, "glacier", "GeoNames $(r[1])", join(filter(!isempty, alt), "|")))
    end
    write_gazetteer("greenland", rows)
end

# SCAR CGA feature types kept, and the label type they map to
const CGA_TYPES = Dict("Glacier" => "glacier", "Ice stream" => "glacier", "Ice shelf" => "iceshelf", "Dome" => "dome")
# the name of a feature is taken from the first of these countries that named it
const CGA_COUNTRIES = ["United States of America", "United Kingdom", "New Zealand", "Australia", "Norway"]

"Antarctic glaciers, ice streams, ice shelves and domes from the SCAR Composite Gazetteer."
function gazetteer_antarctica(g)
    d, h = readdlm(raw("scar", "SCAR_CGA_place_names.csv"), ',', Any; header=true, quotes=true)
    col(n) = d[:, findfirst(==(n), vec(h))]
    name, country, lat, lon = string.(col("place_name_mapping")), string.(col("country_name")), col("latitude"), col("longitude")
    ftype, id = string.(col("feature_type_name")), col("scar_common_id")
    groups = Dict{Any, Vector{Int}}()
    deleted = "is_deleted" in h ? col("is_deleted") .== "Y" : falses(length(name))
    for k in eachindex(name)
        haskey(CGA_TYPES, ftype[k]) && lat[k] isa Real && !deleted[k] && push!(get!(groups, id[k], Int[]), k)
    end
    rows = []
    for ks in values(groups)
        rank(k) = something(findfirst(==(country[k]), CGA_COUNTRIES), length(CGA_COUNTRIES) + 1)
        k = ks[argmin(rank.(ks))]
        inside(g, PROJ_ANT, lon[k], lat[k]) || continue
        alt = filter(!=(name[k]), unique(name[ks]))
        push!(rows, (name[k], Float64(lat[k]), Float64(lon[k]), CGA_TYPES[ftype[k]], "SCAR CGA $(id[k]) ($(country[k]))",
                     join(alt, "|")))
    end
    write_gazetteer("antarctica", rows)
end

rgi_shapefiles() = [joinpath(d, f) for d in readdir(raw("rgi60"); join=true) if isdir(d)
                    for f in readdir(d) if endswith(f, ".shp")]

imbie_shapefile(reg) = only([joinpath(d, f) for d in readdir(raw("imbie"); join=true) if isdir(d) && startswith(basename(d), reg)
                             for f in readdir(d) if endswith(f, ".shp")])

function main(regions)
    for reg in regions
        reg == "greenland"  && prepare_greenland()
        reg == "antarctica" && prepare_antarctica()
    end
    rm(TMP; recursive=true, force=true)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(isempty(ARGS) ? ["greenland", "antarctica"] : ARGS)
end
