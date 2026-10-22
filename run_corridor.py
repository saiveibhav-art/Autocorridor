"""Usage (inside QGIS Python env):  python run_corridor.py  x1 y1  x2 y2  [tolerance_pct]
Coordinates are in the friction raster's CRS (metres)."""
import json, sys
from corridor import run_corridor

cfg = json.load(open("config.json"))
x1, y1, x2, y2 = map(float, sys.argv[1:5])
tol = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
res = run_corridor(cfg["output_friction"], (x1, y1), (x2, y2), "output/corridor", tol)
print(json.dumps(res, indent=2))
