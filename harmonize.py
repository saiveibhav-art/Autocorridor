"""Step 1: put every input on one CRS / extent / pixel grid.

Run inside the QGIS Python console, or standalone after QgsApplication init
and processing.core.Processing.initialize().
"""
import os
import processing
from qgis.core import QgsRasterLayer

# Processing enum for gdal:warpreproject RESAMPLING: 0 nearest, 1 bilinear
NEAREST, BILINEAR = 0, 1


def _warp(src, dst, crs, res, resampling, nodata=-9999, extent=None, extent_crs=None):
    params = {
        "INPUT": src, "TARGET_CRS": crs, "RESAMPLING": resampling,
        "NODATA": nodata, "TARGET_RESOLUTION": res, "DATA_TYPE": 0,
        "OUTPUT": dst,
    }
    if extent:
        params["TARGET_EXTENT"] = extent
        params["TARGET_EXTENT_CRS"] = extent_crs or crs
    return processing.run("gdal:warpreproject", params)["OUTPUT"]


def _rasterize_proximity(vec, ref, crs, res, extent, work, name):
    """Reproject vector, burn to the reference grid, return distance-to-feature raster (m)."""
    proj = processing.run("native:reprojectlayer", {
        "INPUT": vec, "TARGET_CRS": crs, "OUTPUT": "TEMPORARY_OUTPUT"})["OUTPUT"]
    burned = os.path.join(work, f"{name}_burn.tif")
    processing.run("gdal:rasterize", {
        "INPUT": proj, "BURN": 1, "USE_Z": False, "UNITS": 1,
        "WIDTH": res, "HEIGHT": res, "EXTENT": extent, "NODATA": 0,
        "INIT": 0, "DATA_TYPE": 0, "EXTRA": "-at", "OUTPUT": burned})
    dist = os.path.join(work, f"{name}_dist.tif")
    processing.run("gdal:proximity", {
        "INPUT": burned, "BAND": 1, "VALUES": "1", "UNITS": 0,
        "MAX_DISTANCE": 0, "NODATA": -9999, "DATA_TYPE": 5, "OUTPUT": dist})
    return dist


def harmonize(cfg):
    work = cfg["workdir"]
    os.makedirs(work, exist_ok=True)
    crs, res, inp = cfg["target_crs"], cfg["resolution_m"], cfg["inputs"]

    # The DEM defines the reference grid (extent + alignment).
    dem = _warp(inp["dem"], os.path.join(work, "dem.tif"), crs, res, BILINEAR)
    ref = QgsRasterLayer(dem)
    e = ref.extent()
    extent = f"{e.xMinimum()},{e.xMaximum()},{e.yMinimum()},{e.yMaximum()} [{crs}]"

    out = {"dem": dem}
    out["landcover"] = _warp(inp["landcover"], os.path.join(work, "landcover.tif"),
                             crs, res, NEAREST, nodata=0, extent=extent)  # categorical
    out["canopy"] = _warp(inp["canopy"], os.path.join(work, "canopy.tif"),
                          crs, res, BILINEAR, extent=extent)

    out["slope"] = processing.run("gdal:slope", {
        "INPUT": dem, "BAND": 1, "SCALE": 1, "AS_PERCENT": False,
        "COMPUTE_EDGES": True, "OUTPUT": os.path.join(work, "slope.tif")})["OUTPUT"]

    out["road_dist"] = _rasterize_proximity(inp["roads"], dem, crs, res, extent, work, "roads")
    out["settlement_dist"] = _rasterize_proximity(inp["settlements"], dem, crs, res, extent, work, "settle")
    return out
