# Step 0: download all raw data into $CRYOMAPS_DATA/raw (see scripts/paths.jl).
#
# NASA Earthdata files (BedMachine, MEaSUREs) need a ~/.netrc entry:
#   machine urs.earthdata.nasa.gov login <user> password <password>
# Everything else is public. Existing files are skipped; partial downloads resume.
#
# Usage: julia --project=. scripts/fetch_data.jl
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
]

# SCAR Composite Gazetteer of Antarctica (all names, CSV via the AADC web feature service)
const CGA_URL = "https://data.aad.gov.au/geoserver/ows?service=wfs&version=2.0.0&request=GetFeature" *
                "&typeNames=aadc:SCAR_CGA_PLACE_NAMES&outputFormat=csv" *
                "&propertyName=place_name_mapping,country_name,latitude,longitude,feature_type_name,scar_common_id"
const CGA_OUT = joinpath(RAW_DIR, "scar", "SCAR_CGA_place_names.csv")

# ETOPO 2022 30 arc-second surface elevation, regional subset via OPeNDAP
# (fills the Greenland poster beyond the BedMachine domain)
const ETOPO_URL = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/30s/30s_surface_elev_netcdf/" *
                  "ETOPO_2022_v1_30s_N90W180_surface.nc"
const ETOPO_OUT = joinpath(RAW_DIR, "etopo", "ETOPO2022_30s_greenland_subset.nc")
const ETOPO_LATS = (54.0, 88.0)
const ETOPO_LONS = (-115.0, 30.0)

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

function fetch_etopo()
    isfile(ETOPO_OUT) && (println("have ", basename(ETOPO_OUT)); return)
    mkpath(dirname(ETOPO_OUT))
    NCDataset(ETOPO_URL) do ds
        lat = ds["lat"][:]; lon = ds["lon"][:]
        ilat = findfirst(>=(ETOPO_LATS[1]), lat):findlast(<=(ETOPO_LATS[2]), lat)
        ilon = findfirst(>=(ETOPO_LONS[1]), lon):findlast(<=(ETOPO_LONS[2]), lon)
        println("fetching ETOPO ", length(ilon), " x ", length(ilat), " ...")
        z = Matrix{Float32}(undef, length(ilon), length(ilat))
        for c in Iterators.partition(eachindex(ilat), 400)     # chunked to keep requests small
            z[:, c] = coalesce.(ds["z"][ilon, ilat[c]], NaN32)
        end
        NCDataset(ETOPO_OUT*".part", "c") do out
            defVar(out, "lon", lon[ilon], ("lon",); attrib=["units" => "degrees_east", "standard_name" => "longitude"])
            defVar(out, "lat", lat[ilat], ("lat",); attrib=["units" => "degrees_north", "standard_name" => "latitude"])
            defVar(out, "z", z, ("lon", "lat"); deflatelevel=4, attrib=["units" => "m", "_FillValue" => NaN32])
            out.attrib["source"] = ETOPO_URL
            out.attrib["reference"] = "NOAA NCEI (2022): ETOPO 2022 15 Arc-Second Global Relief Model. doi:10.25921/fd45-gt74"
            out.attrib["accessed"] = string(today())
        end
    end
    mv(ETOPO_OUT*".part", ETOPO_OUT)
end

function main()
    for (sub, url, edl) in FILES
        dest = joinpath(RAW_DIR, sub, basename(url))
        download(url, dest; earthdata=edl)
        endswith(dest, ".zip") && unzip(dest)
    end
    download(CGA_URL, CGA_OUT)
    fetch_etopo()
    println("raw data in ", RAW_DIR)
end

main()
