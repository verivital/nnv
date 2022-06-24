import json
import numpy as np

dim = 2
volumes = []

fname1 = f"logged/stat_0037.json"
with open(fname1, "r") as f:
    log1 = json.load(f)

radius1 = log1["stats"]["radius"]
semiAxes_prod1 = log1["stats"]["semiAxes_prod"]
volumes1 = log1["stats"]["volume"]

fname2 = f"logged/stat_0038.json"
with open(fname2, "r") as f:
    log2 = json.load(f)

radius2 = log2["stats"]["radius"]
semiAxes_prod2 = log2["stats"]["semiAxes_prod"]
volumes2 = log2["stats"]["volume"]

for r1, s1, v1, r2, s2, v2 in zip(radius1, semiAxes_prod1, volumes1, radius2, semiAxes_prod2, volumes2):
    d = 2 * max(r1*s1, r2*s2)
    v = d ** dim
    print(v)
    v = min(v, v1)
    v = min(v, v2)
    volumes.append(v)

volumes = np.array(volumes)

print("average volume", volumes.mean())