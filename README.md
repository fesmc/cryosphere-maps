# cryosphere-maps

A0 poster maps of the Greenland (portrait) and Antarctic (landscape) ice sheets:
surface ice velocity over hillshaded topography, drainage regions, surface
contours, the median sea-ice edge, key numbers and place names. Julia +
CairoMakie, with GDAL (via `GDAL_jll`) for regridding. All input data are
downloaded by the pipeline itself (step 0). An interactive version is on
<https://fesmc.github.io/cryosphere-maps/> (step 3). Maps of the ice sheets at
the Last Glacial Maximum are made the same way (see [Paleo maps](#paleo-maps)),
and old-world style maps with expedition routes and wildlife from the same data
(see [Vintage maps](#vintage-maps)).

**Quick start:** the repo contains the poster grids at every 4th node
(`data/prepared/`), enough for the 1600 px sharing PNGs. After the setup below,
`julia --project=. scripts/greenland.jl` makes `plots/greenland_A0_velocity_coarse_small.png`
in ~20 s, without downloading anything. Steps 0–1 are only needed for the
full-resolution posters (`full`) and web tiles.

## Setup

```bash
git clone --single-branch git@github.com:fesmc/cryosphere-maps.git
```
`--single-branch` skips the `gh-pages` branch, which holds the published
website (map tiles and poster PDFs, ~160 MB) and is only needed to publish it.
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

### Full-resolution posters on your own machine

The PDFs and full PNGs are not in the repo (download them from the website,
or make them yourself; they stay untracked). The `pdf` option makes them from
the coarse grids in the repo, without any downloads. Steps 0–1 run on any machine with
Julia, internet access and an Earthdata login; GDAL comes with the Julia
packages. Plotting at full resolution needs ~16 GB of memory.

```bash
export CRYOMAPS_DATA=~/cryomaps-data
```
```bash
julia --project=. scripts/fetch_data.jl
```
```bash
julia --project=. scripts/prepare.jl
```
```bash
julia --project=. scripts/greenland.jl full
```

The SLURM jobs in `jobs/` run the same scripts on albedo.

## Data sources

| Field | Greenland | Antarctica |
|---|---|---|
| Bed/surface topography, thickness, ice mask | BedMachine Greenland v6, 150 m (NSIDC IDBMG4)* | BedMachine Antarctica v4, 500 m (NSIDC-0756)* |
| Surface velocity | MEaSUREs multi-year mosaic, 250 m (NSIDC-0670)* | MEaSUREs InSAR v2, 450 m (NSIDC-0484)* |
| Drainage regions | IMBIE 2 (Rignot & Mouginot): NO, NE, NW, CW, SW, SE | IMBIE 2: 18 basins A-Ap … K-A |
| Beyond BedMachine | ETOPO 2022 30″ (NOAA NCEI, OPeNDAP subset); RGI 6.0 outlines, regions 03/04/06 (OGGM mirror) | — |
| Sea-ice edge | Sea Ice Index v4 (NSIDC G02135), median 1981–2010, March and September | same, September and February |
| Place names | `data/labels_greenland.csv` (curated, GeoNames) | `data/labels_antarctica.csv` (curated, SCAR Composite Gazetteer) |
| Gazetteer (all glacier names) | GeoNames `GL.zip`, feature code GLCR | SCAR Composite Gazetteer (AADC web feature service): glaciers, ice streams, ice shelves, domes |

\* needs an Earthdata login.

The paleo maps also use PaleoMIST 1.0 (Gowan et al., 2021;
doi:[10.1594/PANGAEA.905800](https://doi.org/10.1594/PANGAEA.905800), 3.5 GB)
and, for the Northern Hemisphere, ETOPO 2022 60″ surface and bedrock elevation
north of 24°N (OPeNDAP subsets).

The label CSVs have columns `name, lat, lon, type, source, tier`. Entries
marked `unverified` in `source` are approximate. `tier` 1 = shown on the
posters, 2 = extra detail. The gazetteers (`data/prepared/gazetteer_*.csv`,
tier 3) also list other names of each feature (`alt`). `labels_greenland.csv`
also has Greenlandic names (`kl`, from GeoNames), used by the vintage map.

Overlays (options `routes` and `fauna`, and the vintage maps):

| Overlay | Source |
|---|---|
| Expedition routes | `data/routes_{greenland,antarctica}.csv`, transcribed from the expedition accounts and published positions, with the source of each point (`source`, `quote`, `note`). `track`: `camps` (positions along the route) or `waypoints` (documented places only, drawn as arcs). Rows with `seq` 0 are places off the track |
| Emperor penguin colonies | BAS 2023 colony locations (Fretwell, 2024, UK PDC, doi:10.5285/fb0547e4-d2c1-4580-8c98-182f1da7d9ae) |
| Adélie penguin colonies | MAPPPD (Humphries et al., 2017, www.penguinmap.com), grouped within 100 km |
| Greenland wildlife | Areas Important to Wildlife (GINR/DCE): polar bear denning, musk-ox calving, walrus haul-outs, narwhal summer areas, one point per area, grouped within 100 km |

The wildlife points (`data/prepared/fauna_*.csv`) are made with
`julia --project=. scripts/fetch_data.jl fauna` and
`julia --project=. scripts/prepare_overlays.jl`.

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
   with ETOPO 2022, and RGI marks the ice caps. Takes ~3 min.
   ```bash
   sbatch jobs/prepare.sh
   ```
   → `$CRYOMAPS_DATA/prepared/{greenland_500m,antarctica_1000m}.nc`, with the
   key numbers (area, volume, maximum thickness, sea-level equivalent; for
   Antarctica also East, West and the Peninsula) as global attributes, and in the repo `data/prepared/`: the same grids at
   every 4th node, the gazetteers and the sea-ice edges (GeoJSON). Commit
   these when they change.

2. **Plot**: the full-resolution posters take ~6 min.
   ```bash
   sbatch jobs/plot.sh
   ```
   → `plots/{greenland,antarctica}_A0_velocity.pdf`, which is true A0 size
   (1 unit = 1 pt), plus `.png` (a 100 dpi preview) and `_small.png` (1600 px
   on the long edge, for sharing). Only the small PNGs are tracked in git; the
   PDFs and full PNGs are published on the website (step 3).

   The plot scripts can also be run directly, on any machine:
   ```bash
   julia --project=. scripts/antarctica.jl cmap=ember
   ```
   Without `full`, they only make the small PNG, from the coarse grids in
   `data/prepared/` (~20 s), with `_coarse` in its name.

   **Options:**
   - **Raster style:** `velocity` (default), `surface` or `bed`.
   - **`dark`:** a dark ocean instead of the default light ocean, where deep
     water fades to white.
   - **`nocontours`:** leave out the 500 m surface contours on grounded ice.
   - **`cmap=<name>`:** the velocity colour map:
     - `classic` (default): beige → green → blue → purple → magenta.
     - `ember`: cold to hot, with fast ice glowing orange.
     - `batlow` and `lajolla`: Crameri's perceptually uniform maps.
   - **`tier=2`:** also show the tier-2 names of the label CSV. Names that do
     not fit without overlaps are left out and listed.
   - **`add=<name>;<name>…`:** add names from the label CSV or the gazetteer,
     matched on any of their names, ignoring case and accents. Unknown names
     are reported with the closest matches, and gazetteer names close to a
     label of the CSV are flagged as possible duplicates.
   - **`list`** or **`list=<text>`:** print the available names (containing
     `<text>`) and exit, e.g. `julia --project=. scripts/greenland.jl list=isbrae`.
   - **`routes`**, **`fauna`:** add the expedition routes and the wildlife
     (see the overlays above).
   - **`full`:** full-resolution grid from `$CRYOMAPS_DATA`, and PDF + PNG
     output.
   - **`pdf`:** PDF + PNG output from the coarse grid in the repo (2 km / 4 km
     grid: fine on screen, softer than `full` when printed at A0).

   The default poster goes to `plots/`; any other option combination goes to
   `plots/variants/`, with the options in the file name (for example
   `antarctica_A0_velocity_ember_small.png`). Output from the coarse grid
   has `_coarse` in its name and is not tracked, so it never replaces the
   tracked small PNGs, which come from `full` runs: `jobs/plot.sh` for the
   default posters and `jobs/variants.sh` for the variants in the website
   gallery. Each run reports any overlapping labels.

3. **Website**: a Quarto site (`site/`) with the poster gallery and
   interactive maps (OpenLayers, in the polar stereographic projections),
   published on GitHub Pages.
   ```bash
   sbatch jobs/web.sh
   ```
   → `site/assets/<region>/`: tile pyramids of four map styles from the
   full-resolution grids, GeoJSON layers and the map configuration. Copy them,
   and the poster PDFs and PNGs, to the machine with Quarto, update the
   gallery from `plots/`, then preview or publish:
   ```bash
   rsync -a albedo0.dmawi.de:models/cryosphere-maps/site/assets/ site/assets/
   ```
   ```bash
   rsync -a 'albedo0.dmawi.de:models/cryosphere-maps/plots/*_A0_velocity.*' plots/
   ```
   ```bash
   julia --project=. scripts/web.jl gallery
   ```
   ```bash
   cd site && quarto preview
   ```
   ```bash
   site/publish.sh
   ```
   `site/publish.sh` renders the site and replaces the `gh-pages` branch by a
   single commit, so tiles and PDFs do not accumulate in the repository.
   Without albedo, `julia --project=. scripts/web.jl` builds all assets from
   the coarse grids in the repo, for a local preview.

To change the resolution or map extent, edit `GRIDS` in `scripts/paths.jl`
and rerun steps 1–2.

## Paleo maps

A0 landscape maps of the Northern Hemisphere ice sheets and of Antarctica at a
past time slice, so far the Last Glacial Maximum (20 ka), in the surface and bed
styles of the present-day posters. Only the small PNGs are made
(`plots/paleo/`, tracked); `pdf` gives the PDF too.

**Method** (`scripts/prepare_paleo.jl`): the changes since the time slice in
the PaleoMIST 1.0 reconstruction are applied to present-day topography, ETOPO
2022 for the Northern Hemisphere and BedMachine Antarctica v4 (the repo copy of
the poster grid) for Antarctica:

- bed: `zb(t) = zb_pd - deform(t)`, with `deform` the relative sea-level
  change of PaleoMIST (Earth deformation and sea-level change; from the global
  0.25° grid outside the regional domains);
- surface: `S(t) = max(S_pd, 0) + G[S_pm(t) - S_pm(0)]`, with `S_pm` the
  PaleoMIST surface (paleotopography, 0 over the ocean) and `G` a 15 km
  Gaussian smoothing, which removes the noise of the 5 km PaleoMIST grids;
- ice thickness `S - zb` where PaleoMIST has ice at `t`, floating where it is
  thinner than flotation.

Present-day ice keeps the detail of the base data and ice beyond it follows the
PaleoMIST surface. PaleoMIST reconstructs grounded ice only, so today's ice
shelves beyond its ice are left out. The PaleoMIST regional grids (5 km, North
America incl. Greenland and Iceland, Eurasia, Antarctica) are used in the
projections given in their `projection_info.sh`. The regions of the
Northern Hemisphere numbers (North America, Greenland, Iceland, Eurasia) are
approximate (`nh_regions`).

**Steps:** step 0 downloads PaleoMIST and the ETOPO subsets with the other
data. Then, on albedo (grids and PNGs, ~10 min):
```bash
sbatch jobs/paleo.sh
```
or step by step, on any machine with the raw data:
```bash
julia --project=. scripts/prepare_paleo.jl time=20 nh antarctica
```
→ `data/prepared/paleo_{nh,antarctica}_20ka.nc` in the repo (commit them), with
the surface, bed and mask at the time slice, the present-day mask and the key
numbers (prefix `pd_` for the present). The maps need only these files:
```bash
julia --project=. scripts/paleo.jl nh time=20 surface
```
→ `plots/paleo/nh_20ka_A0_surface_small.png`. Options: `nh` or `antarctica`,
`time=<ka>` (a PaleoMIST time slice, every 2.5 ka; 20 by default), `surface`
(default) or `bed`, and `dark`, `nocontours` and `pdf` as for the present-day
posters. Labels come from `data/labels_paleo_nh.csv` (ice sheets, seas) and,
for Antarctica, the regions and seas of `data/labels_antarctica.csv`.

## Vintage maps

The same data as the posters, drawn as an old chart: inks on paper, ice flow in
rust to oxblood, water lining along the coasts, waves on the sea beyond the
median winter sea-ice edge and floes inside it, expedition routes, wildlife,
a compass rose and a neatline. The Greenland map uses Greenlandic names where
the label CSV has them (`kl`). Fonts: IM Fell English (SIL Open Font License,
`data/fonts/`).
```bash
julia --project=. scripts/vintage.jl antarctica
```
→ `plots/vintage/antarctica_A0_vintage_coarse_small.png`. Options: `greenland`
or `antarctica`, `noroutes`, `nofauna`, and `full` and `pdf` as for the
posters.

## Scripts

- `scripts/paths.jl`: data locations and poster grid definitions.
- `scripts/fetch_data.jl`: step 0, the downloads.
- `scripts/prepare.jl`: step 1, regridding with gdalwarp, ogr2ogr and
  gdal_rasterize, key numbers, gazetteers and sea-ice edges.
- `scripts/greenland.jl`, `scripts/antarctica.jl`: step 2, the poster
  layouts.
- `scripts/common.jl`: hillshade, colour compositing, contours, map
  furniture, command-line options and saving the output.
- `scripts/projection.jl`: polar stereographic projection, its inverse and
  scale factor.
- `scripts/labels.jl`: label styles, name lookup and layouts. Both posters use
  `:coastal`: names are pushed along the coast normal, then slid along the
  coast to avoid overlaps. An alternative `:columns` layout (stacked margin
  columns, as in Rignot & Mouginot 2012) is also available.
- `scripts/prepare_paleo.jl`, `scripts/paleo.jl`, `jobs/paleo.sh`: the paleo
  grids and maps (see [Paleo maps](#paleo-maps)).
- `scripts/overlays.jl`: expedition routes and wildlife glyphs;
  `scripts/prepare_overlays.jl`: the wildlife points.
- `scripts/vintage.jl`: the vintage maps (see [Vintage maps](#vintage-maps)).
- `scripts/smooth.jl`: Gaussian smoothing of grids.
- `scripts/web.jl`: step 3, the web map assets and the poster galleries.
- `site/`: the Quarto website; `site/js/map.js` is the interactive map,
  `site/publish.sh` publishes it.
- `scripts/make_qr.jl`: a one-off that wrote the QR code of the website
  (`data/qr_site.txt`).
- `scripts/add_tiers.jl`: a one-off that added the `tier` column to the label
  CSVs. It has already been applied.

## Known limitations

- **Ice caps outside BedMachine (Greenland):** these come from RGI at 500 m
  and have no velocity data, so they are shown in plain hillshade, without a
  velocity colour.
- **Velocity gaps:** MEaSUREs has gaps in places (for example parts of the
  peripheral glaciers and some ice shelves); these are shown in plain
  hillshade too.
