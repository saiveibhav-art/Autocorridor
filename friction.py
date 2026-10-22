"""Step 2: reclassify each harmonized layer to 1-100 resistance and combine them,
then ADD a crossing penalty on road cells.

Puyravaud et al. (2016) give road cells the highest cost in their elephant model,
so roads are treated as a penalty on the road cells themselves (not a distance decay).
The penalty is added after the weighted average, so friction can exceed 100 on roads.
"""
import csv
import os
import numpy as np
from osgeo import gdal

NODATA = -9999.0


def _read(path):
    ds = gdal.Open(path)
    arr = ds.GetRasterBand(1).ReadAsArray().astype("float32")
    nd = ds.GetRasterBand(1).GetNoDataValue()
    mask = np.zeros(arr.shape, bool) if nd is None else (arr == nd)
    return arr, mask, ds


def _scale(x, lo=1.0, hi=100.0):
    return np.clip(x, lo, hi)


def landcover_cost(arr, lookup_csv, default=50.0):
    table = {int(r["class"]): float(r["resistance"]) for r in csv.DictReader(open(lookup_csv))}
    out = np.full(arr.shape, default, "float32")
    for k, v in table.items():
        out[arr == k] = v
    return _scale(out)


def canopy_cost(pct):
    return _scale(100.0 - np.clip(pct, 0, 100))          # more canopy -> less friction


def slope_cost(deg, flat, steep):
    return _scale(1 + 99 * np.clip((deg - flat) / (steep - flat), 0, 1))


def proximity_cost(dist_m, decay_m):
    return _scale(100.0 * np.exp(-dist_m / decay_m))      # close to feature -> high friction


def road_cells(layers):
    """Boolean raster, True where a road was burned onto the grid (value 1)."""
    burn = os.path.join(os.path.dirname(layers["road_dist"]), "roads_burn.tif")
    arr, _, _ = _read(burn)
    return arr == 1


def build_friction(layers, cfg):
    w = cfg["weights"]
    assert abs(sum(w.values()) - 1.0) < 1e-6, "weights must sum to 1"

    lc, m1, ref = _read(layers["landcover"])
    can, m2, _ = _read(layers["canopy"])
    slp, m3, _ = _read(layers["slope"])
    st, m5, _ = _read(layers["settlement_dist"])

    parts = {
        "landcover": landcover_cost(lc, cfg["landcover_lookup"]),
        "canopy": canopy_cost(can),
        "slope": slope_cost(slp, cfg["slope"]["flat_deg"], cfg["slope"]["max_deg"]),
        "settlement_dist": proximity_cost(st, cfg["decay_m"]["settlement"]),
    }
    friction = sum(w[k] * parts[k] for k in parts).astype("float32")
    friction = friction + float(cfg["road_penalty"]) * road_cells(layers)
    friction[m1 | m2 | m3 | m5] = NODATA

    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(cfg["output_friction"], ref.RasterXSize, ref.RasterYSize, 1,
                    gdal.GDT_Float32, ["COMPRESS=DEFLATE", "TILED=YES"])
    ds.SetGeoTransform(ref.GetGeoTransform())
    ds.SetProjection(ref.GetProjection())
    b = ds.GetRasterBand(1)
    b.SetNoDataValue(NODATA)
    b.WriteArray(friction)
    ds.FlushCache()
    return cfg["output_friction"]
