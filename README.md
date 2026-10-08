# cryosphere-maps

A0 poster maps of the Greenland (portrait) and Antarctic (landscape) ice sheets:
surface ice velocity over hillshaded topography, drainage regions and place
names. Julia + CairoMakie, with GDAL (via `GDAL_jll`) for regridding. All input
data are downloaded by the pipeline itself (step 0).

## Setup

```bash
git clone git@github.com:fesmc/cryosphere-maps.git
```
```bash
cd cryosphere-maps && julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Raw downloads and prepared grids go to `$CRYOMAPS_DATA`. The default is
`/albedo/work/projects/p_forclima/cryosphere-maps-data`; set the variable to
use another location (see `scripts/paths.jl`). About 12 GB are needed.

The NASA Earthdata files need an account and a `~/.netrc` entry (file mode 600):

```
machine urs.earthdata.nasa.gov login <user> password <password>
```

## Data sources

| Field | Greenland | Antarctica |
|---|---|---|
| Bed/surface topography, ice mask | BedMachine Greenland v6, 150 m (NSIDC IDBMG4)* | BedMachine Antarctica v4, 500 m (NSIDC-0756)* |
| Surface velocity | MEaSUREs multi-year mosaic, 250 m (NSIDC-0670)* | MEaSUREs InSAR v2, 450 m (NSIDC-0484)* |
| Drainage regions | IMBIE 2 (Rignot & Mouginot): NO, NE, NW, CW, SW, SE | IMBIE 2: 18 basins A-Ap … K-A |
| Beyond BedMachine | ETOPO 2022 30″ (NOAA NCEI, OPeNDAP subset); RGI 6.0 outlines, regions 03/04/06 (OGGM mirror) | — |
| Place names | `data/labels_greenland.csv` (GeoNames) | `data/labels_antarctica.csv` (SCAR Composite Gazetteer) |

\* needs an Earthdata login.

The label CSVs have columns `name, lat, lon, type, source, tier`. Entries
marked `unverified` in `source` are approximate. `tier` 1 = shown on the
posters, 2 = extra detail.

## Steps

On albedo, run step 0 on the login node (it needs internet access); steps 1–2
run as SLURM jobs from the repo root, with logs in `logs/`.

0. **Download** the raw data to `$CRYOMAPS_DATA/raw`. Existing files are
   skipped and partial downloads resume. This takes ~10 min.
   ```bash
   julia --project=. scripts/fetch_data.jl
   ```

1. **Prepare**: regrid onto the poster grids defined in `scripts/paths.jl`:
   Greenland at 500 m (EPSG:3413), Antarctica at 1 km (EPSG:3031), and the
   basins on a coarser 4 / 8 km grid. Continuous fields are area-averaged and
   the mask uses the mode. Outside the BedMachine domain, Greenland is filled
   with ETOPO 2022, and RGI marks the ice caps. Takes ~2 min.
   ```bash
   sbatch jobs/prepare.sh
   ```
   → `$CRYOMAPS_DATA/prepared/{greenland_500m,antarctica_1000m}.nc`

2. **Plot**: takes ~6 min.
   ```bash
   sbatch jobs/plot.sh
   ```
   → `plots/{greenland,antarctica}_A0_velocity.pdf`, which is true A0 size
   (1 unit = 1 pt), plus `.png` (a 100 dpi preview) and `_small.png` (1600 px
   on the long edge, for sharing).

   **Options** (passed to `jobs/plot.sh` or to the plot scripts):
   - **Raster style:** `velocity` (default), `surface` or `bed`.
   - **`dark`:** a dark ocean instead of the default light ocean, where deep
     water fades to white.
   - **`nocontours`:** leave out the 500 m surface contours on grounded ice.
   - **`cmap=<name>`:** the velocity colour map:
     - `classic` (default): beige → green → blue → purple → magenta.
     - `ember`: cold to hot, with fast ice glowing orange.
     - `batlow` and `lajolla`: Crameri's perceptually uniform maps.
   - **`draft`:** uses every 4th grid node and writes only `*_draft.png`. It
     takes ~1 min, so it can run on the login node while prototyping:
     ```bash
     julia --project=. scripts/antarctica.jl draft cmap=ember
     ```

   Non-default options are appended to the output file name, for example
   `antarctica_A0_velocity_ember.pdf`. Each run reports any overlapping labels.

To change the resolution or map extent, edit `GRIDS` in `scripts/paths.jl`
and rerun steps 1–2.

## Scripts

- `scripts/paths.jl`: data locations and poster grid definitions.
- `scripts/fetch_data.jl`: step 0, the downloads.
- `scripts/prepare.jl`: step 1, regridding with gdalwarp, ogr2ogr and
  gdal_rasterize.
- `scripts/greenland.jl`, `scripts/antarctica.jl`: step 2, the poster
  layouts.
- `scripts/common.jl`: projection, hillshade, colour compositing, contours,
  map furniture, and saving A0 output.
- `scripts/labels.jl`: label styles and layouts. Both posters use `:coastal`:
  names are pushed along the coast normal, then slid along the coast to avoid
  overlaps. An alternative `:columns` layout (stacked margin columns, as in
  Rignot & Mouginot 2012) is also available.
- `scripts/add_tiers.jl`: a one-off that added the `tier` column to the label
  CSVs. It has already been applied.

## Known limitations

- **Ice caps outside BedMachine (Greenland):** these come from RGI at 500 m
  and have no velocity data, so they are shown in plain hillshade, without a
  velocity colour.
- **Velocity gaps:** MEaSUREs has gaps in places (for example parts of the
  peripheral glaciers and some ice shelves); these are shown in plain
  hillshade too.
