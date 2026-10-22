import numpy as np, processing
from osgeo import gdal

lc = r"C:\autocorridor\data\landcover.tif"
binary = r"C:\autocorridor\data\built_binary.tif"
out = r"C:\autocorridor\data\settlements.gpkg"

ds = gdal.Open(lc)
built = (ds.GetRasterBand(1).ReadAsArray() == 50).astype("uint8")

o = gdal.GetDriverByName("GTiff").Create(binary, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Byte)
o.SetGeoTransform(ds.GetGeoTransform())
o.SetProjection(ds.GetProjection())
b = o.GetRasterBand(1)
b.WriteArray(built)
b.SetNoDataValue(0)
o.FlushCache()
o = None

processing.run("gdal:polygonize", {"INPUT": binary, "BAND": 1, "FIELD": "DN",
                                   "EIGHT_CONNECTEDNESS": False, "EXTRA": "",
                                   "OUTPUT": out})
print("done:", out)