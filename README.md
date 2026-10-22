# AutoCorridor

A small PyQGIS/Python tool that turns raw raster and vector layers (elevation, land cover, canopy cover, roads, settlements) into an ecological **resistance (friction) surface** and computes a **least-cost path and corridor swath** between two points. It also includes a one-at-a-time **sensitivity analysis** and a **consensus map** that show how much the result depends on the settings.

> **Status: methods demonstration, not validated.** The resistance values are partly taken from a published elephant study and partly my own choices. Nothing here has been compared with animal movement data. **Do not use the output for conservation or planning decisions.** The point of the project is the workflow and an honest look at how fragile such a model is.

## What it does

1. **Harmonize** (`harmonize.py`): reproject everything to one metric CRS and align all layers to the DEM grid (30 m). Derives slope and a road raster / distance-to-settlement raster. Land cover is resampled with nearest-neighbour (categorical), everything else bilinear.
2. **Friction surface** (`friction.py`): each layer is converted to a 1–100 resistance score, combined as a weighted average, then a **crossing penalty is added on road cells**.
3. **Least-cost path and corridor** (`corridor.py`): Dijkstra on an 8-connected grid; moving between two cells costs the mean of their friction × step length (1 or √2). The corridor swath is every cell whose best start→cell→end route costs at most `(1 + tolerance)` × the least-cost path. Pure NumPy + `heapq`, so no packages need installing inside QGIS.
4. **Sensitivity analysis** (`sensitivity.py`) and **consensus map** (`consensus.py`): re-run the corridor under 16 alternative settings and compare.

## Study area and data (demonstration run)

A ~20 × 20 km box near Bannerghatta National Park, Karnataka, India (77.51–77.69 °E, 12.71–12.89 °N), UTM 43N (EPSG:32643), 30 m cells. One start/end pair, about 12 km apart (UTM 43N: start 790062, 1411027; end 781615, 1420007).

| Layer | Source |
|---|---|
| Elevation | Copernicus GLO-30 DEM (30 m), downloaded via OpenTopography |
| Land cover | ESA WorldCover 2021 v200 (10 m) |
| Roads | OpenStreetMap via QuickOSM, `highway` = motorway to tertiary (+ link roads); 1,561 line features in the box |
| Settlements | Built-up class (50) of WorldCover, polygonized |
| Canopy cover | Hansen Global Forest Change (GFC-2024-v1.12) `treecover2000`, **excluded from the final model** (see below) |

### Why canopy cover was excluded
The Hansen `treecover2000` layer averaged **3.0 %** over cells that WorldCover labels as tree cover in this box, while the same file read **50.6 %** in a test box in the Western Ghats (Kodagu). The file is intact; the layer simply does not describe this area well (it records year-2000 canopy and tall trees only). Its weight was set to 0 instead of letting a near-constant layer add noise.

## Resistance values and what is cited

Land-cover resistances follow **Table 1, "original" column, of Puyravaud et al. (2016)**, an elephant connectivity study in the Nilgiri Biosphere Reserve. The mapping from WorldCover classes to the paper's classes is **my judgement**:

| WorldCover class | Value | Paper class | Confidence |
|---|---|---|---|
| 10 Tree cover | 1 | Evergreen / deciduous forest | clean match |
| 20 Shrubland | 10 | Sparse scrubs and grass | good |
| 30 Grassland | 10 | Sparse scrubs and grass | good |
| 40 Cropland | 30 | Non-irrigated fields (30; irrigated is 50) | **my choice**, WorldCover does not separate irrigated land |
| 50 Built-up | 20 | Habitations / barren land | **likely understates dense urban land** (the paper's area is rural) |
| 60 Bare / sparse | 20 | Barren fields | good |
| 80 Permanent water | 50 | Water bodies | good (paper's alternative model uses 100) |
| 90 Herbaceous wetland | 10 | none | **my assumption**, 0.06 % of the box |

**Not cited (my choices):** the layer weights (land cover 0.606, slope 0.182, settlement distance 0.212, canopy 0), the slope thresholds (5°–30°), the settlement distance decay (1500 m), and the road crossing penalty (100, added after averaging). The paper *adds* its layers' resistances; this tool takes a weighted average, so the result is **informed by** the paper, not a replication of it.

## Results for the demonstration box

**Baseline:** least-cost path 13.51 km (straight line ≈ 12.3 km); 5 % corridor 29.8 km².

### Sensitivity (one change at a time)
`path overlap` = share of the variant's path within 150 m of the baseline path. `corridor Jaccard` = overlap of the two 5 % corridors (1 = identical).

| Variant | Path km | Path overlap % | Corridor km² | Corridor Jaccard | Road cells on path |
|---|---|---|---|---|---|
| **Baseline** | 13.51 | 100 | 29.8 | 1.00 | 8 |
| Road penalty 0 | 13.18 | 86.5 | 23.1 | 0.76 | 18 |
| Road penalty 25 | 13.24 | 86.5 | 25.1 | 0.83 | 11 |
| Road penalty 50 | 13.25 | 86.3 | 27.2 | 0.91 | 10 |
| Road penalty 200 | 13.58 | 51.0 | 29.9 | 0.89 | 6 |
| Road penalty 400 | 13.92 | 39.9 | 32.4 | 0.80 | 5 |
| Land cover weight 0.4 | 13.11 | 84.5 | 29.9 | 0.94 | 7 |
| Land cover weight 0.8 | 14.56 | 44.7 | 30.7 | 0.79 | 6 |
| Built-up 50 | 13.76 | 48.1 | 22.6 | 0.74 | 6 |
| Built-up 95 | 13.85 | 47.9 | 23.2 | 0.71 | 6 |
| Cropland 10 | 13.19 | 85.0 | 28.6 | 0.91 | 7 |
| Cropland 50 | 13.54 | 81.3 | 35.4 | 0.80 | 7 |
| Water 100 | 13.51 | 100 | 29.5 | 0.99 | 8 |
| Settlement decay 750 m | 13.51 | 100 | 29.0 | 0.95 | 8 |
| Settlement decay 3000 m | 13.48 | 50.1 | 30.1 | 0.97 | 6 |
| Slope max 20° | 13.51 | 100 | 29.7 | 0.99 | 8 |
| Slope max 40° | 13.51 | 100 | 29.8 | 1.00 | 8 |

**Reading it:**
- **The single least-cost line is not robust.** 6 of the 16 variants share only 40–51 % of the baseline path. The model is choosing between two routes about 2 km apart (estimated by eye from a map); settings that make it avoid human land more strongly tip it from one to the other.
- **The corridor is steadier than the line.** Corridor overlap stays between 0.71 and 1.00 across variants, because each 5 % corridor is wide enough to contain both routes.
- **Slope barely matters in this box**, and the road penalty cannot remove all crossings (5 road cells remain even at a penalty of 400).

### Consensus map (baseline + 16 variants)
Share of the corridor footprint by number of variants that include a cell:

| Agreement | Area |
|---|---|
| In every variant | 17.8 km² |
| In ≥ 75 % of variants | 26.5 km² |
| In ≥ 50 % of variants | 29.3 km² |
| In at least one variant | 39.2 km² |

This measures robustness to the settings I chose to test, **not** the probability that an animal uses a cell. Cells near the start and end points are in every variant by construction.

### Corridor width depends on the tolerance
In an earlier run (placeholder resistance values, before the cited ones were added), 2 / 5 / 10 / 20 % tolerance gave 12,791 / 25,472 / 47,382 / 120,041 cells (≈ 11.5 / 22.9 / 42.6 / 108 km²). The tolerance has no biological meaning; report corridors as a range, and tie the width to species-specific guidance before using it.

## Limitations

- **No validation** against GPS collars, camera traps, conflict records, or expert-drawn corridors.
- **Resistance values are relative and partly invented** (see above). The class mapping is the weakest link.
- **Least-cost assumes the animal knows the best route.** Circuit-theory methods (Circuitscape, Omniscape) spread movement over many paths and often predict movement better.
- **Point-to-point only**, so results depend on where the two points are placed. Connecting whole habitat patches is the natural next version.
- **Mixed data ages and scales:** WorldCover is 2021 at 10 m resampled to 30 m; narrow roads and small features can vanish or widen. Roads are rasterized "all touched", which makes them wider than they are, and a diagonal step might occasionally slip between road cells.
- **Pure-Python pathfinding** is fine for ~20 km × 20 km at 30 m (a few seconds per path) but would be too slow for state-sized grids.
- **The demonstration box covers the northern, largely urban edge of the Bannerghatta landscape**, not the main elephant corridors.
- **Tested in one environment only:** the standalone QGIS 3.34 install on Windows, one study box, one start/end pair.

## How to run

Requirements: QGIS 3.34 (tested on Windows, standalone install). No extra Python packages are needed.

1. Put your five layers in `data/` as `dem.tif`, `landcover.tif`, `canopy.tif`, `roads.gpkg`, `settlements.gpkg`, and edit `config.json` (CRS, resolution, weights, road penalty).
2. Edit `lookup_landcover.csv` if you use a different land-cover scheme or species.
3. In the QGIS Python Console:

```python
import os, sys
os.chdir(r"C:\autocorridor"); sys.path.insert(0, r"C:\autocorridor")
sys.argv = ["run.py"]
exec(open("run.py").read())                       # builds output/friction_roadpen.tif
```
Remove the friction layer from the QGIS project before rebuilding (Windows locks open files).

4. Corridor between two points (UTM metres, last value = tolerance in %):

```python
from corridor import run_corridor
r = run_corridor("output/friction_roadpen.tif", (790062, 1411027), (781615, 1420007), "output/corridor_t5", 5)
print(r)
```
5. Sensitivity table and consensus map (edit `START_XY` / `END_XY` at the top of `sensitivity.py` for your own points):

```python
exec(open("sensitivity.py", encoding="utf-8").read())
exec(open("consensus.py", encoding="utf-8").read())
```

## Files

| File | Role |
|---|---|
| `harmonize.py` | reproject / align layers, derive slope and road / settlement rasters |
| `friction.py` | reclassify, weight and combine layers; add road crossing penalty |
| `corridor.py` | least-cost path and corridor swath (NumPy + `heapq`) |
| `run.py`, `run_corridor.py` | entry points for the two main steps |
| `sensitivity.py`, `consensus.py` | sensitivity table and consensus map |
| `config.json`, `lookup_landcover.csv` | all tunable settings |

## Possible next steps
1. Replace the remaining uncited settings with sourced, species-specific values.
2. Validate against movement or occurrence data and report a quantitative measure.
3. Connect habitat patches instead of two points; add a circuit-theory comparison.
4. Wrap as a QGIS Processing algorithm.

## References and data attribution
- Puyravaud, J.-P., Cushman, S.A., Davidar, P. & Madappa, D. (2016). Predicting landscape connectivity for the Asian elephant in its largest remaining subpopulation. *Animal Conservation* 20(3): 225–234. doi:10.1111/acv.12314
- ESA WorldCover 2021 v200 (Zanaga et al., 2022), doi:10.5281/zenodo.7254221. Please keep the ESA attribution if you reuse the map.
- Copernicus GLO-30 DEM (ESA / European Union), obtained via OpenTopography.
- Hansen, M.C. et al. (2013). High-resolution global maps of 21st-century forest cover change. *Science* 342: 850–853 (GFC-2024-v1.12 update).
- © OpenStreetMap contributors (ODbL), roads obtained with QuickOSM.

## Author, licence and acknowledgements
- Author: Sai Veibhav
- Licence: MIT
- Developed with assistance from Claude (Anthropic) for code and analysis design; all data checks and decisions are the author's.
