"""Consensus map: re-run the corridor for the baseline and every sensitivity variant, then count,
for each cell, how many variants put it inside the 5% corridor (and how many paths pass near it).

Cells counted in (nearly) every variant are robust to the settings I tested; cells in only a few
are sensitive. This measures robustness to the tested settings, NOT the probability that animals use a cell.

Run from C:\\autocorridor in the QGIS Python Console:
    exec(open("consensus.py", encoding="utf-8").read())
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.getcwd())
import sensitivity as s


def run_all(config_path="config.json"):
    cfg = json.load(open(config_path))
    layers, gt = s.load_layers(cfg)
    roads = layers["roads"][0]
    base = {"lc": s.read_lookup(cfg["landcover_lookup"]), "weights": cfg["weights"],
            "slope": cfg["slope"], "decay": cfg["decay_m"], "road_penalty": cfg["road_penalty"]}
    variants = [("BASELINE", base)] + s.make_variants(base)
    corridor_count = path_count = None
    for i, (name, p) in enumerate(variants, 1):
        pm, sw, ln, res, rc = s.solve(s.build(layers, p), gt, roads)
        thick = s.dilate(pm, 1)                      # 3-cell-wide path, easier to see
        corridor_count = sw.astype("int16") if corridor_count is None else corridor_count + sw
        path_count = thick.astype("int16") if path_count is None else path_count + thick
        print("%2d/%d  %s" % (i, len(variants), name))
    return cfg, gt, corridor_count, path_count, len(variants), res


def write_tif(path, arr, cfg):
    from osgeo import gdal
    ref = gdal.Open(os.path.join(cfg["workdir"], "landcover.tif"))
    t = gdal.GetDriverByName("GTiff").Create(path, ref.RasterXSize, ref.RasterYSize, 1, gdal.GDT_Byte)
    t.SetGeoTransform(ref.GetGeoTransform())
    t.SetProjection(ref.GetProjection())
    b = t.GetRasterBand(1)
    b.WriteArray(arr.astype("uint8"))
    b.SetNoDataValue(0)
    t.FlushCache()


def summarize(corridor_count, n, res):
    km2 = res * res / 1e6
    out = []
    for label, frac in (("in every variant", 1.0), ("in >= 75% of variants", 0.75),
                        ("in >= 50% of variants", 0.5), ("in at least one variant", 0.0)):
        cells = (corridor_count >= max(1, int(np.ceil(frac * n)))).sum()
        out.append((label, int(cells), round(cells * km2, 1)))
    return out


def main():
    cfg, gt, cc, pc, n, res = run_all()
    os.makedirs("output", exist_ok=True)
    write_tif("output/consensus_corridor.tif", cc, cfg)
    write_tif("output/consensus_paths.tif", pc, cfg)
    print("\nCorridor area (5%% tolerance) by agreement across %d variants:" % n)
    for label, cells, a in summarize(cc, n, res):
        print("  %-28s %8d cells  %6.1f km2" % (label, cells, a))
    print("\nSaved: output/consensus_corridor.tif  (cell value = number of variants whose corridor includes it)")
    print("Saved: output/consensus_paths.tif     (cell value = number of variants whose path passes within 1 cell)")


if __name__ == "__main__" or "iface" in globals():
    main()