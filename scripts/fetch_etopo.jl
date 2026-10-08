# Step 1: fetch a regional subset of ETOPO 2022 (60 arc-second surface elevation)
# from NOAA NCEI's THREDDS/OPeNDAP server. Only the subset is transferred.
#
# Output: data/external/ETOPO2022_60s_greenland_subset.nc  (lon, lat, z [m])
# Usage:  julia --project=. scripts/fetch_etopo.jl

using NCDatasets, Dates

const URL = "https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_surface_elev_netcdf/" *
            "ETOPO_2022_v1_60s_N90W180_surface.nc"
const OUT = joinpath(@__DIR__, "..", "data", "external", "ETOPO2022_60s_greenland_subset.nc")

# Region covering the extended Greenland poster window (with some slack)
const LATS = (54.0, 88.0)
const LONS = (-115.0, 30.0)
const STRIDE = 2                      # every 2nd point: 2 arc-min (~3.7 km)

function main()
    mkpath(dirname(OUT))
    NCDataset(URL) do ds
        lat = ds["lat"][:]; lon = ds["lon"][:]
        ilat = findfirst(>=(LATS[1]), lat):STRIDE:findlast(<=(LATS[2]), lat)
        ilon = findfirst(>=(LONS[1]), lon):STRIDE:findlast(<=(LONS[2]), lon)
        println("fetching ", length(ilon), " x ", length(ilat), " points ...")
        z = Float32.(coalesce.(ds["z"][ilon, ilat], NaN32))
        NCDataset(OUT, "c") do out
            defVar(out, "lon", lon[ilon], ("lon",); attrib=["units" => "degrees_east"])
            defVar(out, "lat", lat[ilat], ("lat",); attrib=["units" => "degrees_north"])
            defVar(out, "z", z, ("lon", "lat"); attrib=["units" => "m", "long_name" => "surface elevation"])
            out.attrib["source"] = URL
            out.attrib["subset"] = "lat $(LATS), lon $(LONS), stride $(STRIDE)"
            out.attrib["reference"] = "NOAA National Centers for Environmental Information (2022): " *
                "ETOPO 2022 15 Arc-Second Global Relief Model. doi:10.25921/fd45-gt74"
            out.attrib["accessed"] = string(today())
        end
    end
    println("wrote ", normpath(OUT), " (", round(filesize(OUT)/1e6, digits=1), " MB)")
end

main()
