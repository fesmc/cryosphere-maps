# Step 1: regrid the raw data (step 0) onto the poster grids defined in
# scripts/paths.jl, using GDAL (gdalwarp / ogr2ogr / gdal_rasterize from GDAL_jll).
#
# Output: $CRYOMAPS_DATA/prepared/<region>_<res>m.nc with, on the poster grid,
#   z_srf, z_bed [m], u [m/yr], mask (0 ocean, 1 ice-free land, 2 grounded ice, 3 floating ice)
# and, on the coarser basin grid (xb, yb), basin (integer id, 0 = none; names in
# the `names` attribute).
#
# Usage: julia --project=. scripts/prepare.jl [greenland] [antarctica]
# Takes ~2 min and <3 GB on albedo (jobs/prepare.sh).

using NCDatasets, GDAL_jll, PROJ_jll
include("paths.jl")

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

"Rasterise polygons of vector file(s) `src` onto grid `g`, burning value `burn`."
function rasterize(src, g; dx=g.dx, burn=1, name="rasterized")
    mkpath(TMP)
    reproj = joinpath(TMP, "$(name).gpkg")
    out = joinpath(TMP, "$(name).nc")
    srcs = src isa AbstractVector ? src : [src]
    rm(reproj; force=true)
    for s in srcs
        run(`$(ogr2ogr_exe()) -q -t_srs EPSG:$(g.epsg) -append -nln layer -nlt PROMOTE_TO_MULTI $reproj $s`)
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

function write_prepared(region, g, fields, (basin, basin_names), attrs)
    x, y = axes_km(g); xb, yb = axes_km(g; dx=g.dxb)
    out = prepared_file(region)
    mkpath(dirname(out)); rm(out; force=true)
    NCDataset(out, "c") do ds
        defVar(ds, "x", x, ("x",); attrib=["units" => "km"])
        defVar(ds, "y", y, ("y",); attrib=["units" => "km"])
        defVar(ds, "xb", xb, ("xb",); attrib=["units" => "km"])
        defVar(ds, "yb", yb, ("yb",); attrib=["units" => "km"])
        for (k, (v, units)) in fields
            T = k == "mask" ? Int8 : Float32
            defVar(ds, k, T.(replace(v, NaN => (T == Int8 ? -1 : NaN32))), ("x", "y"); deflatelevel=4,
                   attrib=["units" => units])
        end
        defVar(ds, "basin", Int16.(basin), ("xb", "yb"); deflatelevel=4,
               attrib=["units" => "1", "names" => join(basin_names, ",")])
        ds.attrib["grid"] = "EPSG:$(g.epsg), $(g.dx) km (basins $(g.dxb) km)"
        for (k, v) in attrs
            ds.attrib[k] = v
        end
    end
    println("wrote ", out)
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
    m  = warp(bm_var(bm, "mask"), g; s_srs=s, resample="mode", name="grl_mask")
    vx = warp(raw("measures", "greenland_vel_mosaic250_vx_v1.tif"), g; s_srs=s, name="grl_vx")
    vy = warp(raw("measures", "greenland_vel_mosaic250_vy_v1.tif"), g; s_srs=s, name="grl_vy")
    et = warp(bm_var(raw("etopo", "ETOPO2022_30s_greenland_subset.nc"), "z"), g;
              s_srs="EPSG:4326", resample="bilinear", name="grl_etopo")
    rgi = rasterize(rgi_shapefiles(), g; name="grl_rgi")

    out = isnan.(m)                                   # outside the BedMachine domain
    zb[out] .= et[out]; zs[out] .= max.(et[out], 0.0)
    m[out] .= ifelse.(et[out] .< 0, 0.0, 1.0)
    m[(m .== 1) .& (rgi .> 0)] .= 2                   # ice caps outside the ice sheet

    basin = rasterize_basins(imbie_shapefile("GRE"), g, "SUBREGION1"; exclude=["ICE_CAP"], name="grl_basin")
    write_prepared("greenland", g,
        ["z_srf" => (zs, "m"), "z_bed" => (zb, "m"), "u" => (speed(vx, vy), "m/yr"), "mask" => (m, "1")], basin,
        ["sources" => "BedMachine Greenland v6; MEaSUREs NSIDC-0670 v1; ETOPO 2022 30s; RGI 6.0; IMBIE2 basins"])
end

"Antarctica: BedMachine v4, MEaSUREs velocity v2 (NSIDC-0484), IMBIE2 basins."
function prepare_antarctica()
    g  = GRIDS["antarctica"]; s = "EPSG:3031"
    bm = raw("bedmachine", "NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc")
    zs = warp(bm_var(bm, "surface"), g; s_srs=s, name="ant_surface")
    zb = warp(bm_var(bm, "bed"), g; s_srs=s, name="ant_bed")
    m  = warp(bm_var(bm, "mask"), g; s_srs=s, resample="mode", name="ant_mask")
    m[m .== 4] .= 2                                   # Lake Vostok -> grounded ice
    m[isnan.(m)] .= 0
    vel = raw("measures", "antarctica_ice_velocity_450m_v2.nc")
    vx = warp(bm_var(vel, "VX"), g; s_srs=s, name="ant_vx")
    vy = warp(bm_var(vel, "VY"), g; s_srs=s, name="ant_vy")

    basin = rasterize_basins(imbie_shapefile("ANT"), g, "Subregion"; name="ant_basin")
    write_prepared("antarctica", g,
        ["z_srf" => (zs, "m"), "z_bed" => (zb, "m"), "u" => (speed(vx, vy), "m/yr"), "mask" => (m, "1")], basin,
        ["sources" => "BedMachine Antarctica v4; MEaSUREs NSIDC-0484 v2; IMBIE2 basins"])
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

main(isempty(ARGS) ? ["greenland", "antarctica"] : ARGS)
