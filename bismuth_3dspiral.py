import Rhino.Geometry as rg
import random
import math
import json

# ---------------------------------------------------------------------------------
# Inputs (keep your existing curve-generation inputs if you want)
# Required minimal inputs for this version:
# start_point, seed, layers, seg_count, initial_len, scale, z_step_cells,
# branch_interval, branch_prob
#
# New inputs for geometry build:
# terrace_w (float), terrace_h (float), terrace_steps (int), shrink (float)
#
# Outputs: B (Breps), C (Curves), L (Lines), N (int)
# ---------------------------------------------------------------------------------


# HELPER FUNCTIONS
def rot90_xy(d): return (-d[1], d[0], 0)
def add(a,b): return (a[0]+b[0], a[1]+b[1], a[2]+b[2])
def mul(d,k): return (d[0]*k, d[1]*k, d[2]*k)

# INPUTS


#params#single GH input: a JSON string in a Panel
if not params:
    params = "{}"

try:
    PDICT = json.loads(params)
except Exception:
    #params#if JSON is broken, show error nicely
    raise ValueError("params#JSON parse failed. Check commas/quotes in your Panel text.")

def P(key, default):
    v = PDICT.get(key, default)
    return default if v is None else v

#params#read values (only here, no GH inputs spam)
start_point = P("start_point", [0,0,0])
start_point = rg.Point3d(start_point[0], start_point[1], start_point[2])

seed = int(P("seed", 1))
layers = int(P("layers", 10))
seg_count = int(P("seg_count", 80))
initial_len = float(P("initial_len", 20.0))
scale = float(P("scale", 0.95))
z_step_cells = float(P("z_step_cells", 3.0))
branch_interval = int(P("branch_interval", 12))
branch_prob = float(P("branch_prob", 0.25))

terrace_w = float(P("terrace_w", 2.0))
terrace_h = float(P("terrace_h", 1.0))
terrace_steps = int(P("terrace_steps", 8))
shrink = float(P("shrink", 0.92))

branch_steps_mul = float(P("branch_steps_mul", 0.35))
branch_shrink_mul = float(P("branch_shrink_mul", 1.06))
branch_kill_p = float(P("branch_kill_p", 0.15))
branch_kill_after = int(P("branch_kill_after", 3))


"""
#old
if start_point is None: start_point = rg.Point3d(0,0,0)
seed = int(seed) if seed is not None else 1
layers = int(layers) if layers else 10
seg_count = int(seg_count) if seg_count else 80
initial_len = float(initial_len) if initial_len else 20.0
scale = float(scale) if scale is not None else 0.95
z_step_cells = float(z_step_cells) if z_step_cells else 3.0
branch_interval = int(branch_interval) if branch_interval else 12
branch_prob = float(branch_prob) if branch_prob is not None else 0.25

terrace_w = float(terrace_w) if terrace_w is not None else 2.0
terrace_h = float(terrace_h) if terrace_h is not None else 1.0
terrace_steps = int(terrace_steps) if terrace_steps else 8
shrink = float(shrink) if shrink is not None else 0.92
"""


random.seed(seed)


#these are the 4 vectors we use later (we might want to add some for ARTISTIC IRREGULARITIES)
# 0 x, 1 y, 2 -x, 3 -y
# right now we are keeping these to lock the direction to 90 degree turns
DIRS2D = [rg.Vector3d(1,0,0), rg.Vector3d(0,1,0), rg.Vector3d(-1,0,0), rg.Vector3d(0,-1,0)]


# HERE WE DRAW A SQUARE SPIRAL POLYLINE
# THE LENGTH * SCALE LINE IS VERY IMPORTANT
# FOR OTHER VARIATIONS THAT SCALE VALUE SHOULD BE A FUNCTION

#from tuning:
#scale closer to 1.0 → spiral tightens slowly (bigger spirals)
#scale smaller (0.90–0.94) → tightens fast (more nested look)

def create_spiral_polyline(pt0, dir0, seg_count, length0, scale, z):
    pts = []
    #current point (start) with original x and y but z is forced
    cur = rg.Point3d(pt0.X, pt0.Y, z)
    #next point vector
    d = rg.Vector3d(dir0)
    #length is the var that gets shorter every turn
    length = float(length0)
    #start point added before loop
    pts.append(cur)

    for i in range(seg_count):
        #makes sure L is at least 0.5 at all times
        L = max(0.5, length)
        #defines the point nxt by adding d.X/Y (part from the vector d) times the L
        nxt = rg.Point3d(
            cur.X + d.X * L, 
            cur.Y + d.Y * L, 
            z
            )
        #adds them to pts
        pts.append(nxt)
        #KEY STEP: Updating the new point as the start point
        cur = nxt
        #scale comes from tuning input (gh)
        length *= scale
        #rotate 90 deg / if 1/0/0 then 0/-1/0, if 0/-1/0 then -1/0/0, if -1/0/0 then 0/1/0 and then back again
        d = rg.Vector3d(d.Y,-d.X, 0)
        #stopper
        if length < 0.8:
            break

    pl = rg.Polyline(pts)
    return pl


#makes a closed boundary from the curves by offsetting +w -w
def offset_closed_outline(polyline, w):
    """
    Create closed outline curve from polyline by offsetting left/right and joining.
    """
    #Rhino’s Offset() works on a Curve, not on a Polyline object directly (or not reliably).
    #so we turn it into a Rhino Curve
    crv = polyline.ToNurbsCurve()
    #define the offset plane by putting in the first point and the z part of vector 
    #for a normal vector (this would be different in more complex directions with var z)
    plane = rg.Plane(rg.Point3d(0,0,polyline[0].Z), rg.Vector3d.ZAxis)
    #Here we use the Offset function from Rhino
    #0.01 is the tolerance, w is direction (in this case just both sides(depends on crv directio))
    off1 = crv.Offset(plane,  w, 0.01, rg.CurveOffsetCornerStyle.Sharp)
    off2 = crv.Offset(plane, -w, 0.01, rg.CurveOffsetCornerStyle.Sharp)
    
    if not off1 or not off2: 
        return None

    #take first offset results (just a convention with this function)
    #and then flips second crv to create a closes Loop (A->B + B->A = ABA) 
    c1 = off1[0]
    c2 = off2[0]
    c2.Reverse()
    #join...
    joined = rg.Curve.JoinCurves([c1, c2], 0.05)
    if not joined:
        return None
    #again, rhino function stuff, calling first works
    out = joined[0]
    #this is actually practical, makes higher tolerance IF NOT CLOSED...
    if not out.IsClosed:
        out.MakeClosed(0.05)
    return out


def planar_brep_from_closed(crv):
    # THIS TAKES A CLOSED PLANAR CURVE AND MAKES A FLAT SURFACE (BREP) FROM IT
    # closed outline -> "cap" surface
    # We are NOT using it in the final extrusion pipeline right now (we extrude curves directly),
    # but it’s useful if later we want planar plates, booleans, or face operations.

    # Rhino function returns a LIST of breps (could be multiple if curve has weird loops)
    breps = rg.Brep.CreatePlanarBreps(crv, 0.05)

    # if we got something, take the first one (convention)
    if breps and len(breps) > 0:
        return breps[0]
    return None


def scale_curve_2d(crv, s):
    #THIS SCALES A CURVE IN 2D AROUND ITS "CENTER"
    #we compute centroid of the area enclosed by the curve
    #then we scale the curve around that centroid

    # rea mass properties only works properly if curve is planar + closed
    amp = rg.AreaMassProperties.Compute(crv)
    if not amp:
        
        #if curve is not computable, return original curve (safe fallback)
        return crv

    #centroid point (the "center" of the area)
    cen = amp.Centroid

    #build a scaling transform around centroid
    xform = rg.Transform.Scale(cen, s)

    #duplicate curve so we don't modify original
    c = crv.DuplicateCurve()
    c.Transform(xform)
    return c


# C will store ALL spiral curves (main spirals + branch spirals)
# L will store line segments 
#(debug view: shows each step as a straight line)
C = []
L = []

# base start point for the first layer
base = start_point

# length_layer is the "starting step length" per layer
# and it shrinks each layer (so higher layers usually get smaller spirals)
length_layer = initial_len


for layer in range(layers):
    # compute Z height of this layer
    # layer 0: z = base.Z + 0 * z_step_cells
    # layer 1: z = base.Z + 1 * z_step_cells
    z = base.Z + layer * z_step_cells

    # choose direction from DIRS2D:
    # layer % 4 cycles: 0,1,2,3,0,1,2,3 ...
    # this keeps a nice "swirl / rotation" feeling across stacked layers
    dir0 = DIRS2D[layer % 4]  # swirl-preserving

    # create the polyline spiral on that z plane
    pl = create_spiral_polyline(base, dir0, seg_count, length_layer, scale, z)

    # convert polyline to a Rhino curve so later functions (offset etc.) work well
    crv = pl.ToNurbsCurve()

    # store it
    C.append(crv)

    # DEBUG LINES: store each polyline segment as a Line for visualization
    for i in range(pl.Count - 1):
        L.append(rg.Line(pl[i], pl[i + 1]))

    
    # OPTIONAL BRANCH SYSTEM
    
    # idea: along the main spiral, at some indices, start mini-spirals
    # This creates "multi-crystal / cluster" vibes, but too much makes "tower city"
    if branch_prob > 0:
        # iterate over vertices of the polyline (excluding first and last)
        for i in range(1, pl.Count - 1):

            # branch condition:
            # every branch_interval-th point, roll a random chance
            if i % branch_interval == 0 and random.random() < branch_prob:

                # direction of the current local segment
                bdir = rg.Vector3d(pl[i + 1] - pl[i])
                bdir.Unitize()

                # force it to one of the axis directions (so branches remain 90° locked)
                # choose the dominant axis:
                if abs(bdir.X) > abs(bdir.Y):
                    # mostly x direction
                    bdir = rg.Vector3d(1 if bdir.X > 0 else -1, 0, 0)
                else:
                    # mostly y direction
                    bdir = rg.Vector3d(0, 1 if bdir.Y > 0 else -1, 0)

                # create a smaller spiral starting from that branch point
                # seg_count*0.25 = fewer segments
                # length_layer*0.55 = shorter start length than main spiral
                bpl = create_spiral_polyline(
                    pl[i],
                    bdir,
                    int(seg_count * 0.25),
                    length_layer * 0.55,
                    scale,
                    z
                )

                # store the branch curve too
                C.append(bpl.ToNurbsCurve())

                # debug lines for branch
                for j in range(bpl.Count - 1):
                    L.append(rg.Line(bpl[j], bpl[j + 1]))

    # shrink the layer’s initial length as we go up
    # max(3.0, ...) prevents it from becoming too tiny
    length_layer = max(3.0, length_layer * 0.86)


#build terraces from curves
#B will store final Breps (extruded terrace plates)
B = []

for crv in C:
    #GOAL:
    #curve (spiral centerline) -> closed ribbon outline -> inset steps -> extrude each step
    #This is our "materialization" system.

    #Step 1) Convert curve to polyline for stable offsets
    #TryGetPolyline() returns (success, polyline) in Python
    ok, pl = crv.TryGetPolyline()

    if not ok:
        #fallback: manually sample the curve into points
        #(because offsetting a wild nurbs curve can cause weird results)
        t0, t1 = crv.Domain.T0, crv.Domain.T1
        pts = [crv.PointAt(t0)]
        stepsN = 120  #more steps = closer approximation but heavier
        for k in range(1, stepsN + 1):
            t = t0 + (t1 - t0) * (k / float(stepsN))
            pts.append(crv.PointAt(t))
        pl = rg.Polyline(pts)

    #Step 2) create a closed outline ribbon from polyline
    #offset +w and -w then join into a boundary
    outline = offset_closed_outline(pl, terrace_w)
    if outline is None:
        # offset failed -> skip this curve
        continue

    #Step 3) make multiple inset outlines (terrace steps)
    #outline_0 = original ribbon boundary
    #utline_1 = scaled inward by shrink
    #outline_2 = scaled inward by shrink^2 ...
    outlines = []
    for s in range(terrace_steps):
        c = outline if s == 0 else scale_curve_2d(outline, (shrink ** s))
        outlines.append(c)

    #Step 4) extrude each outline and stack them upwards
    #This gives a stepped terrace volume (like bismuth)
    for s, oc in enumerate(outlines):
        if oc is None or not oc.IsClosed:
            continue

        #move outline up by s * terrace_h
        #so each step sits above the previous one
        ocT = oc.DuplicateCurve()
        ocT.Transform(rg.Transform.Translation(0, 0, s * terrace_h))

        #extrude the closed curve into a solid "plate"
        #Extrusion.Create(CURVE, HEIGHT, CAP)
        ext = rg.Extrusion.Create(ocT, terrace_h, True)

        #if extrusion succeeded, convert to Brep and store it
        if ext:
            B.append(ext.ToBrep())

#N = how many breps we ended with
N = len(B)
