import Rhino.Geometry as rg
import math
import random



# --- 2. THE REUSABLE FUNCTION ---
def grow_three_segments(start_pt, start_dir, length, shrink_factor, is_clockwise):
    lines = []
    curr_pos = start_pt
    curr_dir = rg.Vector3d(start_dir)
    curr_len = length
    
    # Toggle rotation: -90 for CW, +90 for CCW
    angle = -math.pi / 2 if is_clockwise else math.pi / 2 
    
    for i in range(9):
        curr_len *= shrink_factor
        curr_dir.Rotate(angle, rg.Vector3d.ZAxis)
        
        new_p = curr_pos + (curr_dir * curr_len)
        lines.append(rg.Line(curr_pos, new_p))
       
        curr_pos = new_p
    #    if distance(curr_pos, new_p) >
    #        angle = -math.pi / 2 if is_clockwise else math.pi / 2
        
    return lines

def run():
        
       import Rhino.Geometry as rg
import random
import math

# ------------------------------------------------------------
# Inputs expected (make GH inputs with these names):
# S (Point or list[Point]), cell (float), steps (int),
# vertical_every (int), branch_p (float),
# min_dist (int), seed (int),
# size_jitter (float), age_penalty (float),
# nutrient_init (float), nutrient_consume (float), nutrient_diffuse (float)
# Outputs: B, P, L
# ------------------------------------------------------------

# ---------- input hygiene ----------
def as_points(x):
    if x is None: return []
    if isinstance(x, list) or isinstance(x, tuple): return list(x)
    return [x]

S = as_points(S)

if not S:
    # default seed
    S = [rg.Point3d(0,0,0)]

cell = float(cell) if cell else 1.0
steps = int(steps) if steps else 3000
vertical_every = max(1, int(vertical_every) if vertical_every else 4)
branch_p = float(branch_p) if branch_p is not None else 0.22
min_dist = int(min_dist) if min_dist is not None else 1
seed = int(seed) if seed is not None else 1
size_jitter = float(size_jitter) if size_jitter is not None else 0.15
age_penalty = float(age_penalty) if age_penalty is not None else 0.015
nutrient_init = float(nutrient_init) if nutrient_init is not None else 1.0
nutrient_consume = float(nutrient_consume) if nutrient_consume is not None else 0.08
nutrient_diffuse = float(nutrient_diffuse) if nutrient_diffuse is not None else 0.06

random.seed(seed)

# ---------- grid utils ----------
def to_grid(pt):
    # snap world point to integer cell coords
    return (int(round(pt.X / cell)), int(round(pt.Y / cell)), int(round(pt.Z / cell)))

def to_world(g):
    return rg.Point3d(g[0]*cell, g[1]*cell, g[2]*cell)

# 6-neighborhood
DIRS = [
    ( 1, 0, 0),
    (-1, 0, 0),
    ( 0, 1, 0),
    ( 0,-1, 0),
    ( 0, 0, 1),
    ( 0, 0,-1),
]

# anisotropic base weights (tune these)
BASE_W = {
    ( 1, 0, 0): 0.36,
    ( 0, 1, 0): 0.26,
    ( 0, 0, 1): 0.18,
    (-1, 0, 0): 0.10,
    ( 0,-1, 0): 0.07,
    ( 0, 0,-1): 0.03,
}

def weighted_choice(items, weights):
    s = sum(weights)
    if s <= 1e-12:
        return None
    r = random.random() * s
    acc = 0.0
    for it, w in zip(items, weights):
        acc += w
        if r <= acc:
            return it
    return items[-1]

# ---------- occupancy + neighbor inhibition ----------
occ = set()                # occupied grid cells
age = {}                   # cell -> age int
nutr = {}                  # cell -> nutrient float
parent = {}                # cell -> parent cell
active_front = set()       # subset of occ; candidates to grow from

def near_occupied(g, r):
    # Chebyshev ball for speed; "good enough" for inhibition
    gx,gy,gz = g
    for dx in range(-r, r+1):
        for dy in range(-r, r+1):
            for dz in range(-r, r+1):
                if (gx+dx, gy+dy, gz+dz) in occ:
                    return True
    return False

# ---------- initialize seeds ----------
for pt in S:
    g0 = to_grid(pt)
    if g0 not in occ:
        occ.add(g0)
        age[g0] = 0
        nutr[g0] = nutrient_init
        parent[g0] = None
        active_front.add(g0)

# ---------- growth scoring ----------
def dir_allowed(d, step_idx):
    # vertical lag: only allow +Z and -Z every N steps
    if d == (0,0,1) or d == (0,0,-1):
        return (step_idx % vertical_every) == 0
    return True

def exposed_faces(g):
    # a cell is a "front" if it has at least one empty neighbor
    gx,gy,gz = g
    for d in DIRS:
        ng = (gx+d[0], gy+d[1], gz+d[2])
        if ng not in occ:
            return True
    return False

def score_dir(src, d, step_idx):
    # Weight components:
    # - base anisotropy
    # - nutrient at source
    # - age penalty (older grows less)
    # - inhibition (avoid crowding)
    w = BASE_W.get(d, 0.0)

    n = nutr.get(src, 0.0)
    a = age.get(src, 0)

    # aging reduces likelihood gradually
    w *= max(0.0, (1.0 - a * age_penalty))

    # nutrient scales growth capacity
    w *= max(0.0, n)

    # vertical lag already handled outside
    # inhibition: penalize growth into crowded space
    gx,gy,gz = src
    tgt = (gx+d[0], gy+d[1], gz+d[2])
    if tgt in occ:
        return 0.0
    if min_dist > 0 and near_occupied(tgt, min_dist):
        w *= 0.12  # strong penalty, but not zero → keeps terraces dense
    return w

# ---------- main loop ----------
links = []  # (src_world, tgt_world)

for step_idx in range(steps):
    if not active_front:
        break

    # pick a source from active fronts, biased by nutrient (more alive grows more)
    fronts = list(active_front)
    fw = [max(1e-6, nutr.get(f, 0.0)) for f in fronts]
    src = weighted_choice(fronts, fw)
    if src is None:
        break

    # if src isn't actually exposed anymore, drop it
    if not exposed_faces(src):
        if src in active_front: active_front.remove(src)
        continue

    # choose direction
    dirs = []
    ws = []
    for d in DIRS:
        if not dir_allowed(d, step_idx):
            continue
        dirs.append(d)
        ws.append(score_dir(src, d, step_idx))

    dsel = weighted_choice(dirs, ws)
    if dsel is None:
        # src is starved or blocked
        nutr[src] = max(0.0, nutr.get(src, 0.0) - nutrient_consume*0.5)
        if nutr[src] <= 1e-3:
            active_front.discard(src)
        continue

    gx,gy,gz = src
    tgt = (gx+dsel[0], gy+dsel[1], gz+dsel[2])

    # place new cell
    if tgt not in occ:
        occ.add(tgt)
        age[tgt] = 0
        parent[tgt] = src
        nutr[tgt] = nutr.get(src, nutrient_init) * 0.85  # inherit but slightly reduced

        # consume nutrient at source (starvation mechanism)
        nutr[src] = max(0.0, nutr.get(src, 0.0) - nutrient_consume)

        # aging update
        age[src] = age.get(src, 0) + 1

        # diffusion: gently share nutrient to neighbors (creates terrace continuity)
        if nutrient_diffuse > 0:
            share = nutr[tgt] * nutrient_diffuse
            nutr[tgt] -= share
            # distribute to any occupied neighbors
            neigh = []
            for d in DIRS:
                ng = (tgt[0]+d[0], tgt[1]+d[1], tgt[2]+d[2])
                if ng in occ:
                    neigh.append(ng)
            if neigh:
                per = share / float(len(neigh))
                for ng in neigh:
                    nutr[ng] = nutr.get(ng, 0.0) + per

        # keep new cell active; source might remain active depending on nutrient
        active_front.add(tgt)
        if nutr[src] <= 1e-3:
            active_front.discard(src)

        links.append((to_world(src), to_world(tgt)))

        # branching: occasionally re-activate src to sprout multiple terraces
        if random.random() < branch_p and nutr.get(src, 0.0) > 0.15:
            active_front.add(src)

# ---------- build geometry ----------
B = []
P = []
L = []

# create jittered boxes for the iridescent "stepped plate" feel
# IMPORTANT: still snapped to grid; jitter only changes box size slightly
for g in occ:
    cpt = to_world(g)
    P.append(cpt)

    j = 1.0 + (random.random()*2.0 - 1.0) * size_jitter
    sx = cell * j
    sy = cell * j
    sz = cell * j

    # box aligned to world axes
    x = rg.Interval(cpt.X - sx*0.5, cpt.X + sx*0.5)
    y = rg.Interval(cpt.Y - sy*0.5, cpt.Y + sy*0.5)
    z = rg.Interval(cpt.Z - sz*0.5, cpt.Z + sz*0.5)
    box = rg.Box(rg.Plane.WorldXY, x, y, z)
    brep = box.ToBrep()
    if brep:
        B.append(brep)

for a,b in links:
    L.append(rg.Line(a,b))
    
    # return outputs
    # return B, P, L


        if __name__ == "__main__":
            t = run(seed=1)
            print(f"terraces: {len(t)}")
            
            
         