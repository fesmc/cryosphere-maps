# Step 0: download all raw data into $CRYOMAPS_DATA/raw (see scripts/paths.jl).
#
# NASA Earthdata files (BedMachine, MEaSUREs) need a ~/.netrc entry:
#   machine urs.earthdata.nasa.gov login <user> password <password>
# Everything else is public. Existing files are skipped; partial downloads resume.
#
# Usage: julia --project=. scripts/fetch_data.jl [fauna]
# `fauna`: only the small wildlife files for the map overlays (scripts/prepare_overlays.jl).
# Run on a node with internet access (albedo login node).

using NCDatasets, Dates
include("paths.jl")

const NSIDC  = "https://data.nsidc.earthdatacloud.nasa.gov/nsidc-cumulus-prod-protected"
const SEAICE = "https://noaadata.apps.nsidc.org/NOAA/G02135"

# (subdirectory, url, needs Earthdata login)
const FILES = [
    # BedMachine Greenland v6 (Morlighem et al.), 150 m, 2.8 GB
    ("bedmachine", "$NSIDC/ICEBRIDGE/IDBMG4/6/1993/01/01/BedMachineGreenland-v6.nc", true),
    # BedMachine Antarctica v4 (Morlighem et al.), 500 m, 1.1 GB
    ("bedmachine", "$NSIDC/MEASURES/NSIDC-0756/4/1970/01/01/NSIDC-0756_BedMachineAntarctica_19700101-20191001_V04.1.nc", true),
    # MEaSUREs multi-year Greenland velocity mosaic (Joughin et al., 2018), 250 m
    ("measures", "$NSIDC/MEASURES/NSIDC-0670/1/1995/12/01/greenland_vel_mosaic250_vx_v1.tif", true),
    ("measures", "$NSIDC/MEASURES/NSIDC-0670/1/1995/12/01/greenland_vel_mosaic250_vy_v1.tif", true),
    # MEaSUREs InSAR-based Antarctica velocity v2 (Rignot et al., 2011/2017), 450 m, 6.5 GB
    ("measures", "$NSIDC/MEASURES/NSIDC-0484/2/1996/01/01/antarctica_ice_velocity_450m_v2.nc", true),
    # IMBIE2 drainage basins (Zwally et al., 2012)
    ("imbie", "http://imbie.org/wp-content/uploads/2016/09/GRE_Basins_IMBIE2_v1.3.zip", false),
    ("imbie", "http://imbie.org/wp-content/uploads/2016/09/ANT_Basins_IMBIE2_v1.6.zip", false),
    # RGI 6.0 outlines (OGGM mirror) for ice caps around Greenland outside BedMachine
    ("rgi60", "https://cluster.klima.uni-bremen.de/~oggm/rgi/www.glims.org/RGI/rgi60_files/03_rgi60_ArcticCanadaNorth.zip", false),
    ("rgi60", "https://cluster.klima.uni-bremen.de/~oggm/rgi/www.glims.org/RGI/rgi60_files/04_rgi60_ArcticCanadaSouth.zip", false),
    ("rgi60", "https://cluster.klima.uni-bremen.de/~oggm/rgi/www.glims.org/RGI/rgi60_files/06_rgi60_Iceland.zip", false),
    # Sea Ice Index v4 (Fetterer et al., NSIDC G02135): median ice edge 1981-2010,
    # March (max) and September (min) in the north, February (min) and September (max) in the south
    [("seaice", "$SEAICE/north/monthly/shapefiles/shp_median/median_extent_N_$(m)_1981-2010_polyline_v4.0.zip", false)
     for m in ("03", "09")]...,
    [("seaice", "$SEAICE/south/monthly/shapefiles/shp_median/median_extent_S_$(m)_1981-2010_polyline_v4.0.zip", false)
     for m in ("02", "09")]...,
    # GeoNames place names for Greenland (glacier names for the gazetteer)
    ("geonames", "https://download.geonames.org/export/dump/GL.zip", false),
    # PaleoMIST 1.0 ice sheet reconstruction, 80-0 ka every 2.5 kyr (Gowan et al., 2021;
    # doi:10.1594/PANGAEA.905800), 3.5 GB, for the paleo maps
    ("paleomist", "https://hs.pangaea.de/Maps/Global_Ice_Sheets/Gowan_ice_reconstruction.zip", false),
]

# SCAR Composite Gazetteer of Antarctica (all names, CSV via the AADC web feature service)
const CGA_URL = "https://data.aad.gov.au/geoserver/ows?service=wfs&version=2.0.0&request=GetFeature" *
                "&typeNames=aadc:SCAR_CGA_PLACE_NAMES&outputFormat=csv" *
                "&propertyName=place_name_mapping,country_name,latitude,longitude,feature_type_name,scar_common_id"
const CGA_OUT = joinpath(RAW_DIR, "scar", "SCAR_CGA_place_names.csv")

# ETOPO 2022 subsets via OPeNDAP (NOAA NCEI, doi:10.25921/fd45-gt74): the 30 arc-second
# surface elevation around Greenland (fills the Greenland poster beyond the BedMachine
# domain) and the 60 arc-second surface and bedrock elevation north of 24°N (present-day
# base of the Northern Hemisphere paleo map).
const ETOPO = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022"
const ETOPO_SUBSETS = [
    ("$ETOPO/30s/30s_surface_elev_netcdf/ETOPO_2022_v1_30s_N90W180_surface.nc", "ETOPO2022_30s_greenland_subset.nc",
     (54.0, 88.0), (-115.0, 30.0)),
    ("$ETOPO/60s/60s_surface_elev_netcdf/ETOPO_2022_v1_60s_N90W180_surface.nc", "ETOPO2022_60s_surface_nh.nc",
     (24.0, 90.0), (-180.0, 180.0)),
    ("$ETOPO/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc", "ETOPO2022_60s_bed_nh.nc",
     (24.0, 90.0), (-180.0, 180.0)),
]

# Wildlife for the map overlays (scripts/prepare_overlays.jl), as (url, file name in raw/fauna):
# emperor penguin colonies 2023 (Fretwell, 2024, UK PDC, doi:10.5285/fb0547e4-d2c1-4580-8c98-182f1da7d9ae),
# MAPPPD penguin counts (Humphries et al., 2017, www.penguinmap.com) and the Greenland Areas
# Important to Wildlife (GINR/DCE): musk-ox calving, walrus haul-outs, polar bear denning and
# narwhal summer areas.
const AIW = "https://services-eu1.arcgis.com/0uK40YtWoUkQMlYW/arcgis/rest/services/Areas_Important_to_Wildlife_data/FeatureServer"
const FAUNA_FILES = [
    ("https://ramadda.data.bas.ac.uk/repository/entry/get/emperor_colony_locations2023.kmz?entryid=" *
     "synth%3Afb0547e4-d2c1-4580-8c98-182f1da7d9ae%3AL2VtcGVyb3JfY29sb255X2xvY2F0aW9uczIwMjMua216",
     "emperor_colony_locations2023.kmz"),
    ("https://www.penguinmap.com/mapppd/DownloadAll/", "mapppd_AllCounts.csv"),
    [("$AIW/$id/query?where=1%3D1&outFields=*&outSR=4326&f=geojson", "aiw_$(name).geojson")
     for (id, name) in ((14, "muskox_calving"), (17, "walrus_haulout"), (21, "polarbear_denning"),
                        (22, "narwhal_summer"))]...,
]

function download(url, dest; earthdata=false)
    if isfile(dest)
        println("have ", basename(dest)); return
    end
    mkpath(dirname(dest))
    part = dest*".part"
    cookies = joinpath(RAW_DIR, ".edl_cookies")
    auth = earthdata ? `-n -c $cookies -b $cookies` : ``
    println("downloading ", basename(dest), " ...")
    run(`curl -sS -L --fail --retry 3 -C - $auth -o $part $url`)
    mv(part, dest)
end

function unzip(zip)
    dir = splitext(zip)[1]
    isdir(dir) && return
    run(`unzip -q -o $zip -d $dir`)
end

"Subset (lats, lons) of the ETOPO 2022 grid at `url`, written to raw/etopo/`name`."
function fetch_etopo(url, name, lats, lons)
    out = joinpath(RAW_DIR, "etopo", name)
    isfile(out) && (println("have ", name); return)
    mkpath(dirname(out))
    NCDataset(url) do ds
        lat = ds["lat"][:]; lon = ds["lon"][:]
        ilat = findfirst(>=(lats[1]), lat):findlast(<=(lats[2]), lat)
        ilon = findfirst(>=(lons[1]), lon):findlast(<=(lons[2]), lon)
        println("fetching ", name, " ", length(ilon), " x ", length(ilat), " ...")
        z = Matrix{Float32}(undef, length(ilon), length(ilat))
        nc = max(1, 8_000_000 ÷ length(ilon))                 # chunked to keep requests small
        for c in Iterators.partition(eachindex(ilat), nc)
            z[:, c] = coalesce.(ds["z"][ilon, ilat[c]], NaN32)
        end
        NCDataset(out*".part", "c") do o
            defVar(o, "lon", lon[ilon], ("lon",); attrib=["units" => "degrees_east", "standard_name" => "longitude"])
            defVar(o, "lat", lat[ilat], ("lat",); attrib=["units" => "degrees_north", "standard_name" => "latitude"])
            defVar(o, "z", z, ("lon", "lat"); deflatelevel=4, attrib=["units" => "m", "_FillValue" => NaN32])
            o.attrib["source"] = url
            o.attrib["reference"] = "NOAA NCEI (2022): ETOPO 2022 15 Arc-Second Global Relief Model. doi:10.25921/fd45-gt74"
            o.attrib["accessed"] = string(today())
        end
    end
    mv(out*".part", out)
end

function main(args=ARGS)
    for (url, name) in FAUNA_FILES
        download(url, joinpath(RAW_DIR, "fauna", name))
    end
    "fauna" in args && return
    for (sub, url, edl) in FILES
        dest = joinpath(RAW_DIR, sub, basename(url))
        download(url, dest; earthdata=edl)
        endswith(dest, ".zip") && unzip(dest)
    end
    download(CGA_URL, CGA_OUT)
    for (url, name, lats, lons) in ETOPO_SUBSETS
        fetch_etopo(url, name, lats, lons)
    end
    println("raw data in ", RAW_DIR)
end

main()
