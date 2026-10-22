"""One-at-a-time sensitivity analysis for AutoCorridor.

Reuses the layers already built by run.py (output/work/*.tif), so it does NOT
re-run harmonize.py. For each variant it rebuilds the friction surface in memory,
finds the least-cost path and the corridor swath, and compares them with the baseline.

Run from C:\\autocorridor inside the QGIS Python Console:
    exec(open("sensitivity.py", encoding="utf-8").read())
"""
import copy
import csv
import json
import math
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.getcwd())
import friction as fr
from corridor import prepare_cost, least_cost_path, corridor_swath, xy_to_rc

START_XY = (790062, 1411027)   # change to your own start/end (UTM metres)
END_XY = (781615, 1420007)
TOL = 5                        # corridor tolerance in percent
NEAR_CELLS = 5                 # "same place" = within 5 cells (150 m at 30 m)


def load_layers(cfg):
    names = {"landcover": "landcover.tif", "canopy": "canopy.tif", "slope": "slope.tif",
             "road_dist": "roads_dist.tif", "settlement_dist": "settle_dist.tif"}
    layers, gt = {}, None
    for key, fname in names.items():
        arr, mask, ds = fr._read(os.path.join(cfg["workdir"], fname))
        layers[key] = (arr, mask)
        gt = ds.GetGeoTransform()
    return layers, gt


def read_lookup(path):
    return {int(r["class"]): float(r["resistance"]) for r in csv.DictReader(open(path))}


def reweight(weights, key, new):
    """Set weights[key]=new and rescale the others so everything still sums to 1."""
    w = dict(weights)
    others = sum(v for k, v in w.items() if k != key)
    for k in w:
        w[k] = new if k == key else w[k] * (1 - new) / others
    return w


def build(layers, p):
    lc, m1 = layers["landcover"]
    can, m2 = layers["canopy"]
    slp, m3 = layers["slope"]
    rd, m4 = layers["road_dist"]
    st, m5 = layers["settlement_dist"]
    lcc = np.full(lc.shape, 50.0, "float32")       # same default as friction.py
    for k, v in p["lc"].items():
        lcc[lc == k] = v
    parts = {
        "landcover": fr._scale(lcc),
        "canopy": fr.canopy_cost(can),
        "slope": fr.slope_cost(slp, p["slope"]["flat_deg"], p["slope"]["max_deg"]),
        "road_dist": fr.proximity_cost(rd, p["decay"]["road"]),
        "settlement_dist": fr.proximity_cost(st, p["decay"]["settlement"]),
    }
    f = sum(p["weights"][k] * parts[k] for k in parts).astype("float32")
    f[m1 | m2 | m3 | m4 | m5] = fr.NODATA
    return f


def dilate(m, r):
    H, W = m.shape
    pad = np.pad(m, r)
    out = np.zeros_like(m)
    for dr in range(-r, r + 1):
        for dc in range(-r, r + 1):
            if dr * dr + dc * dc <= r * r:
                out |= pad[r + dr: r + dr + H, r + dc: r + dc + W]
    return out


def solve(f, gt):
    cost = prepare_cost(f)
    start, end = xy_to_rc(gt, *START_XY), xy_to_rc(gt, *END_XY)
    path, _ = least_cost_path(cost, start, end)
    _, swath = corridor_swath(cost, start, end, TOL)
    pm = np.zeros(cost.shape, bool)
    rr, cc = zip(*path)
    pm[list(rr), list(cc)] = True
    res = abs(gt[1])
    length = sum(math.hypot(b[0] - a[0], b[1] - a[1]) for a, b in zip(path, path[1:])) * res
    return pm, swath, length, res


def make_variants(base):
    def v(name, **change):
        p = copy.deepcopy(base)
        for k, val in change.items():
            if k == "lc":
                p["lc"].update(val)
            elif k == "weights":
                p["weights"] = val
            elif k == "slope":
                p["slope"].update(val)
            elif k == "decay":
                p["decay"].update(val)
        return name, p
    w = base["weights"]
    return [
        v("landcover weight 0.4", weights=reweight(w, "landcover", 0.4)),
        v("landcover weight 0.6", weights=reweight(w, "landcover", 0.6)),
        v("built-up 50", lc={50: 50}),
        v("built-up 95 (old placeholder)", lc={50: 95}),
        v("cropland 10", lc={40: 10}),
        v("cropland 50", lc={40: 50}),
        v("water 100", lc={80: 100}),
        v("road decay 250 m", decay={"road": 250}),
        v("road decay 1000 m", decay={"road": 1000}),
        v("settlement decay 750 m", decay={"settlement": 750}),
        v("settlement decay 3000 m", decay={"settlement": 3000}),
        v("slope max 20 deg", slope={"max_deg": 20}),
        v("slope max 40 deg", slope={"max_deg": 40}),
    ]


def main(config_path="config.json", out_csv="output/sensitivity.csv"):
    cfg = json.load(open(config_path))
    layers, gt = load_layers(cfg)
    base = {"lc": read_lookup(cfg["landcover_lookup"]), "weights": cfg["weights"],
            "slope": cfg["slope"], "decay": cfg["decay_m"]}
    t0 = time.time()
    b_path, b_sw, b_len, res = solve(build(layers, base), gt)
    b_near = dilate(b_path, NEAR_CELLS)
    rows = [["variant", "path_km", "path_within_150m_of_base_pct", "base_within_150m_of_path_pct",
             "corridor_km2", "corridor_overlap_jaccard"]]
    rows.append(["BASELINE", round(b_len / 1000, 2), 100.0, 100.0,
                 round(b_sw.sum() * res * res / 1e6, 1), 1.0])
    print("baseline done (%.0f s)" % (time.time() - t0))
    for name, p in make_variants(base):
        pm, sw, ln, _ = solve(build(layers, p), gt)
        near_v = dilate(pm, NEAR_CELLS)
        rows.append([name, round(ln / 1000, 2),
                     round(100 * (pm & b_near).sum() / pm.sum(), 1),
                     round(100 * (b_path & near_v).sum() / b_path.sum(), 1),
                     round(sw.sum() * res * res / 1e6, 1),
                     round((sw & b_sw).sum() / (sw | b_sw).sum(), 2)])
        print("%-32s done (%.0f s)" % (name, time.time() - t0))
    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="") as f:
        csv.writer(f).writerows(rows)
    print()
    for r in rows:
        print("%-32s %8s %10s %10s %10s %8s" % tuple(r))
    print("\nSaved:", out_csv)
    return rows


if __name__ == "__main__" or "iface" in globals():
    main()