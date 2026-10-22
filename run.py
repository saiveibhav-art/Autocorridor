import json, sys
from harmonize import harmonize
from friction import build_friction

cfg = json.load(open(sys.argv[1] if len(sys.argv) > 1 else "config.json"))
layers = harmonize(cfg)
print("Friction surface written to:", build_friction(layers, cfg))
