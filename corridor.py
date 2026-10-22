"""Step 3: least-cost path + corridor swath on a friction surface.

No extra packages needed: core logic uses only NumPy (bundled with QGIS) and
Python's built-in heapq. File I/O uses GDAL/OGR, which also ship with QGIS.
"""
import heapq
import math
import numpy as np

NODATA = -9999.0
_SQ2 = math.sqrt(2.0)
# 8-connected neighbours: (d_row, d_col, step length in cells)
_NEIGHBOURS = [(-1, 0, 1.0), (1, 0, 1.0), (0, -1, 1.0), (0, 1, 1.0),
               (-1, -1, _SQ2), (-1, 1, _SQ2), (1, -1, _SQ2), (1, 1, _SQ2)]


# ---------- core (no GIS dependencies) ----------
def prepare_cost(friction, nodata=NODATA):
    """Impassable cells (nodata) become inf; everything else must be > 0."""
    cost = friction.astype("float64").copy()
    cost[(cost == nodata) | ~np.isfinite(cost)] = np.inf
    if np.nanmin(cost[np.isfinite(cost)]) <= 0:
        raise ValueError("Friction values must be > 0")
    return cost


def _dijkstra(cost, origin, target=None, track_path=False):
    """Dijkstra over an 8-connected grid. Moving between two cells costs
    mean(friction_a, friction_b) x step length (1 or sqrt2).
    Stops early at `target` if given."""
    rows, cols = cost.shape
    flat = cost.ravel().tolist()
    dist = [math.inf] * (rows * cols)
    prev = [-1] * (rows * cols) if track_path else None
    s = origin[0] * cols + origin[1]
    t = -1 if target is None else target[0] * cols + target[1]
    dist[s] = 0.0
    heap = [(0.0, s)]
    pop, push = heapq.heappop, heapq.heappush
    while heap:
        d, u = pop(heap)
        if d > dist[u]:
            continue
        if u == t:
            break
        r, c = divmod(u, cols)
        cu = flat[u]
        for dr, dc, w in _NEIGHBOURS:
            rr, cc = r + dr, c + dc
            if rr < 0 or rr >= rows or cc < 0 or cc >= cols:
                continue
            v = rr * cols + cc
            cv = flat[v]
            if cv == math.inf:
                continue
            nd = d + 0.5 * (cu + cv) * w
            if nd < dist[v]:
                dist[v] = nd
                if track_path:
                    prev[v] = u
                push(heap, (nd, v))
    return dist, prev


def least_cost_path(cost, start, end):
    """start/end are (row, col). Returns (list of (row, col), total cost)."""
    rows, cols = cost.shape
    dist, prev = _dijkstra(cost, start, end, track_path=True)
    e = end[0] * cols + end[1]
    if dist[e] == math.inf:
        raise ValueError("No route exists between start and end (blocked by NoData)")
    path, u = [], e
    while u != -1:
        path.append(divmod(u, cols))
        u = prev[u]
    return path[::-1], dist[e]


def cost_distance(cost, origin):
    """Accumulated cost from one origin cell to every cell."""
    dist, _ = _dijkstra(cost, origin)
    return np.array(dist).reshape(cost.shape)


def corridor_swath(cost, start, end, tolerance_pct=10.0):
    """Cells whose best route through them (start -> cell -> end) costs at most
    (1 + tolerance_pct/100) x the least-cost path. Returns (sum_cost_surface, mask)."""
    total = cost_distance(cost, start) + cost_distance(cost, end)
    best = np.nanmin(total[np.isfinite(total)])
    mask = total <= best * (1 + tolerance_pct / 100.0)
    return total, mask


# ---------- georeferencing helpers ----------
def xy_to_rc(gt, x, y):
    col = int((x - gt[0]) / gt[1])
    row = int((y - gt[3]) / gt[5])
    return row, col


def rc_to_xy(gt, row, col):
    return gt[0] + (col + 0.5) * gt[1], gt[3] + (row + 0.5) * gt[5]


# ---------- file I/O ----------
def run_corridor(friction_tif, start_xy, end_xy, out_dir, tolerance_pct=10.0):
    import os
    from osgeo import gdal, ogr, osr

    os.makedirs(out_dir, exist_ok=True)
    ds = gdal.Open(friction_tif)
    gt, proj = ds.GetGeoTransform(), ds.GetProjection()
    band = ds.GetRasterBand(1)
    friction = band.ReadAsArray()
    cost = prepare_cost(friction, band.GetNoDataValue() or NODATA)

    start, end = xy_to_rc(gt, *start_xy), xy_to_rc(gt, *end_xy)
    for name, (r, c) in (("start", start), ("end", end)):
        if not (0 <= r < cost.shape[0] and 0 <= c < cost.shape[1]) or not np.isfinite(cost[r, c]):
            raise ValueError(f"{name} point falls outside the raster or on a NoData cell")

    path, total = least_cost_path(cost, start, end)
    swath, mask = corridor_swath(cost, start, end, tolerance_pct)

    # corridor raster: accumulated cost inside the swath, NoData elsewhere
    out = np.where(mask, swath, NODATA).astype("float32")
    tif = os.path.join(out_dir, "corridor.tif")
    o = gdal.GetDriverByName("GTiff").Create(tif, ds.RasterXSize, ds.RasterYSize, 1,
                                             gdal.GDT_Float32, ["COMPRESS=DEFLATE"])
    o.SetGeoTransform(gt); o.SetProjection(proj)
    ob = o.GetRasterBand(1); ob.SetNoDataValue(NODATA); ob.WriteArray(out); o.FlushCache()

    # least-cost path as a line in a GeoPackage
    gpkg = os.path.join(out_dir, "lcp.gpkg")
    drv = ogr.GetDriverByName("GPKG")
    if os.path.exists(gpkg):
        drv.DeleteDataSource(gpkg)
    vds = drv.CreateDataSource(gpkg)
    srs = osr.SpatialReference(); srs.ImportFromWkt(proj)
    lyr = vds.CreateLayer("least_cost_path", srs, ogr.wkbLineString)
    lyr.CreateField(ogr.FieldDefn("cost", ogr.OFTReal))
    lyr.CreateField(ogr.FieldDefn("length_m", ogr.OFTReal))
    line = ogr.Geometry(ogr.wkbLineString)
    for r, c in path:
        line.AddPoint_2D(*rc_to_xy(gt, r, c))
    feat = ogr.Feature(lyr.GetLayerDefn())
    feat.SetField("cost", float(total)); feat.SetField("length_m", line.Length())
    feat.SetGeometry(line); lyr.CreateFeature(feat); vds = None

    return {"corridor_raster": tif, "path_layer": gpkg, "path_cost": float(total),
            "path_length_m": line.Length(), "corridor_cells": int(mask.sum())}
