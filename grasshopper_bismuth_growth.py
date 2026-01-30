"""
Grasshopper Python3 component script
Inputs (suggested):
  C: Center Point (Point3d)
  S: Seed Points (list of Point3d)
  N: Steps (int)
  Step: Step size (float)
  Attract: Attraction strength (float)
  Jitter: Random jitter strength (float)
  Starve: Starve radius (float)
  Seed: Random seed (int)
Outputs:
  Curves: list of Polylines
  Points: list of lists of Point3d
"""

import Rhino.Geometry as rg
import random

# -------------------------
# Defaults / input hygiene
# -------------------------
if C is None:
    C = rg.Point3d(0, 0, 0)

if S is None:
    S = [rg.Point3d(20, 0, 0), rg.Point3d(-20, 0, 0), rg.Point3d(0, 20, 0), rg.Point3d(0, -20, 0)]

if N is None:
    N = 120

if Step is None:
    Step = 1.0

if Attract is None:
    Attract = 1.0

if Jitter is None:
    Jitter = 0.35

if Starve is None:
    Starve = 6.0

if Seed is None:
    Seed = 1

rng = random.Random(Seed)

# -------------------------
# Growth function
# -------------------------

def jitter_vector(scale, local_rng):
    vec = rg.Vector3d(
        local_rng.uniform(-1.0, 1.0),
        local_rng.uniform(-1.0, 1.0),
        local_rng.uniform(-0.35, 0.35),
    )
    if vec.IsZero:
        return rg.Vector3d(0, 0, 0)
    vec.Unitize()
    return vec * scale


def grow_path(seed_pt, index):
    pts = [seed_pt]
    local_rng = random.Random(rng.randint(0, 10**9) + index * 9973)

    for _ in range(int(N)):
        last = pts[-1]
        to_center = C - last
        dist = to_center.Length

        # Starve before reaching the center
        if dist <= Starve:
            break

        if not to_center.IsZero:
            to_center.Unitize()

        # Slightly vary attraction per seed to showcase different expressions
        attract_scale = Attract * (0.75 + 0.5 * (index % 5) / 4.0)

        # Combine attraction and jitter
        direction = to_center * attract_scale + jitter_vector(Jitter, local_rng)

        if direction.IsZero:
            break

        direction.Unitize()
        next_pt = last + direction * Step
        pts.append(next_pt)

    return pts

# -------------------------
# Build outputs
# -------------------------
all_points = []
curves = []

for i, seed in enumerate(S):
    pts = grow_path(seed, i)
    all_points.append(pts)
    if len(pts) > 1:
        curves.append(rg.Polyline(pts))
    else:
        curves.append(rg.Polyline([seed]))

Points = all_points
Curves = curves
