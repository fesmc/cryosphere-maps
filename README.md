# cryosphere-maps

A0 poster maps of the Greenland (portrait) and Antarctic (landscape) ice sheets:
surface ice velocity over hillshaded topography, drainage divides and place
names. Julia + CairoMakie.

## Setup

```bash
git clone git@github.com:fesmc/cryosphere-maps.git
```
```bash
cd cryosphere-maps && julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The ice-sheet fields are read from a local `~/models/ice_data` checkout (path
set by `ICE_DATA` in `scripts/common.jl`).

## Input data

**Local (`~/models/ice_data`)**, 8 km unless noted:

| Field | Greenland | Antarctica |
|---|---|---|
| Topography | `GRL-8KM/GRL-8KM_TOPO-M17.nc` (BedMachine v3) | `ANT-8KM/ANT-8KM_TOPO-BedMachine.nc` (BedMachine v2) |
| Velocity | `GRL-8KM/GRL-8KM_VEL-J18.nc` | `ANT-16KM/ANT-16KM_VEL-R11-2.nc` (16 km) |
| Basins | `GRL-8KM/GRL-8KM_BASINS-nasa.nc` | `ANT-16KM/ANT-16KM_BASINS-nasa.nc` (16 km) |

The Antarctic velocity and basins come from the 16 km files because the 8 km
versions are broken: `ANT-8KM_VEL-R11-2.nc` has `uxy_srf` and `uy_srf` all
zero, and `ANT-8KM_BASINS-nasa.nc` has `basin` all zero.

**Downloaded** (Greenland only; raw files go to `data/external/`, which is not
committed; the regridded results in `data/GRL-8KM-EXT_*.nc` are). These fill the map beyond the BedMachine
domain (x −720…960 km):
- ETOPO 2022 topography/bathymetry (NOAA NCEI, doi:10.25921/fd45-gt74), steps 1–2.
- RGI 6.0 glacier outlines (RGI Consortium, 2017, doi:10.7265/N5-RGI-60) for
  regions 03 Arctic Canada North, 04 Arctic Canada South and 06 Iceland. These
  mask the ice caps; steps 3–4.

**Place names:** `data/labels_{greenland,antarctica}.csv` (`name, lat, lon,
type, source, tier`). They were compiled from the SCAR Composite Gazetteer
(Antarctica), GeoNames (Greenland) and Wikipedia. Entries marked `unverified`
in `source` are approximate. `tier` 1 = shown on the posters, 2 = extra detail.

## Steps

1. **Fetch the ETOPO 2022 subset** (Greenland only). This is a ~18 MB subset
   (54–88°N, 115°W–30°E, 2 arc-min) read over OPeNDAP from NOAA THREDDS; the
   global file is not downloaded.
   ```bash
   julia --project=. scripts/fetch_etopo.jl
   ```
   → `data/external/ETOPO2022_60s_greenland_subset.nc`

2. **Regrid ETOPO onto the extended 8 km Greenland grid.** This is the GRL-8KM
   grid padded by 400 km on each side, filled by bin-averaging the projected
   ETOPO points.
   ```bash
   julia --project=. scripts/prepare_etopo.jl
   ```
   → `data/GRL-8KM-EXT_ETOPO2022.nc`

3. **Fetch the RGI 6.0 outlines.** These are regional shapefiles from the OGGM
   mirror (no login): 03 (11 MB), 04 (26 MB) and 06 (2.5 MB).
   ```bash
   julia --project=. scripts/fetch_rgi.jl
   ```
   → `data/external/rgi60/<region>/*.shp`

4. **Rasterise RGI onto the same extended grid** as a glacier area fraction
   (4×4 sub-samples per 8 km cell). Cells with a fraction > 0.5 are treated as ice.
   ```bash
   julia --project=. scripts/prepare_rgi.jl
   ```
   → `data/GRL-8KM-EXT_RGI60.nc`

5. **Plot**:
   ```bash
   julia --project=. scripts/greenland.jl
   ```
   ```bash
   julia --project=. scripts/antarctica.jl
   ```
   → `plots/{greenland,antarctica}_A0_velocity.{pdf,png}`. The PDF is true
   A0 size (1 unit = 1 pt); the PNG is a 100 dpi preview. Pass `surface` or
   `bed` as an argument for the other raster styles.

## Scripts

- `scripts/common.jl`: projection, regridding/refinement, hillshade,
  colour compositing, contours, map furniture, and saving A0 output.
- `scripts/labels.jl`: label styles and layouts. Both posters use `:coastal`:
  names are pushed along the coast normal, then slid along the coast to avoid
  overlaps. An alternative `:columns` layout (stacked margin columns, as in
  Rignot & Mouginot 2012) is also available.
- `scripts/add_tiers.jl`: a one-off that added the `tier` column to the label
  CSVs. It has already been applied.

## Known limitations

- **Resolution:** the source fields are 8 km (16 km for Antarctic velocity and
  basins). They are bilinearly refined for smooth rendering, but this adds no
  detail.
- **Greenland outside BedMachine:** ice comes from RGI at 8 km (fraction > 0.5),
  so small glaciers are dropped. There is no velocity data there, so these ice
  caps are shown in plain hillshade, without a velocity colour.
