import Rhino.Geometry as rg
import random
import math
import json

# ---------------------------------------------------------------------------------
# Inputs:
#params#single GH input: a JSON string in a Panel
#
# Outputs (wire these as GH outputs):
#B#Breps#terrace solids
#C#Curves#spiral centerlines (main+branches)
#L#Lines#debug segment lines
#N#int#number of breps
#DBG#tuple#failure counters
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

seed = int(P("seed", 9))
layers = int(P("layers", 13))
seg_count = int(P("seg_count", 30))
initial_len = float(P("initial_len", 45.0))
scale = float(P("scale", 0.85))
z_step_cells = float(P("z_step_cells", 2.7))
#Smaller intervalls >>> More Branching
branch_interval = int(P("branch_interval", 7))
#High branching probality vs low
branch_prob = float(P("branch_prob", 0.85))

terrace_w = float(P("terrace_w", 0.8))
terrace_h = float(P("terrace_h", 0.8))
terrace_steps = int(P("terrace_steps", 8))
shrink = float(P("shrink", 0.92))

#params#branch "die out" controls
#first one growth length z
branch_steps_mul = float(P("branch_steps_mul", 1.75))
branch_shrink_mul = float(P("branch_shrink_mul", 1.03))
branch_kill_p = float(P("branch_kill_p", 0.15))
branch_kill_after = int(P("branch_kill_after", 15))

#params#optional area stop to prevent needle pillars (units^2, depends on your model scale)
min_step_area = float(P("min_step_area", 1.5))

random.seed(seed)

#these are the 4 vectors we use later (we might want to add some for ARTISTIC IRREGULARITIES)
#0 x, 1 y, 2 -x, 3 -y
#right now we are keeping these to lock the direction to 90 degree turns
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


def offset_closed_outline(polyline, w):
    #make ribbon outline from open polyline: offset left/right then cap ends and join
    crv = polyline.ToNurbsCurve()
    plane = rg.Plane(rg.Point3d(0,0,polyline[0].Z), rg.Vector3d.ZAxis)

    off1 = crv.Offset(plane,  w, 0.01, rg.CurveOffsetCornerStyle.Sharp)
    off2 = crv.Offset(plane, -w, 0.01, rg.CurveOffsetCornerStyle.Sharp)
    if not off1 or not off2:
        return None

    c1 = off1[0]
    c2 = off2[0]

    #important#make both curves run in same "loop direction"
    c2r = c2.DuplicateCurve()
    c2r.Reverse()

    #cap#connect ends so we ALWAYS get a closed boundary
    a0 = c1.PointAtStart
    a1 = c1.PointAtEnd
    b0 = c2r.PointAtStart
    b1 = c2r.PointAtEnd

    cap1 = rg.LineCurve(a1, b0)#cap at one end
    cap2 = rg.LineCurve(b1, a0)#cap at the other end

    joined = rg.Curve.JoinCurves([c1, cap1, c2r, cap2], 0.05)
    if not joined:
        return None

    out = joined[0]
    if not out.IsClosed:
        out.MakeClosed(0.05)

    return out


def planar_brep_from_closed(crv):
    # THIS TAKES A CLOSED PLANAR CURVE AND MAKES A FLAT SURFACE (BREP) FROM IT
    # closed outline -> "cap" surface
    # We are NOT using it in the final extrusion pipeline right now (we extrude curves directly),
    # but it’s useful if later we want planar plates, booleans, or face operations.
    breps = rg.Brep.CreatePlanarBreps(crv, 0.05)
    if breps and len(breps) > 0:
        return breps[0]
    return None


def scale_curve_2d(crv, s):
    #THIS SCALES A CURVE IN 2D AROUND ITS "CENTER"
    #we compute centroid of the area enclosed by the curve
    #then we scale the curve around that centroid
    amp = rg.AreaMassProperties.Compute(crv)
    if not amp:
        #if curve is not computable, return original curve (safe fallback)
        return crv
    cen = amp.Centroid
    xform = rg.Transform.Scale(cen, s)
    c = crv.DuplicateCurve()
    c.Transform(xform)
    return c


# C will store ALL spiral curves (main spirals + branch spirals)
# L will store line segments
#(debug view: shows each step as a straight line)
C = []
L = []
Cmeta = []#new#parallel list aligned to C so we know main vs branch in terrace stage

# base start point for the first layer
base = start_point

# length_layer is the "starting step length" per layer
# and it shrinks each layer (so higher layers usually get smaller spirals)
length_layer = initial_len


for layer in range(layers):
    # compute Z height of this layer
    z = base.Z + layer * z_step_cells

    # choose direction from DIRS2D:
    dir0 = DIRS2D[layer % 4]  # swirl-preserving

    # create the polyline spiral on that z plane
    pl = create_spiral_polyline(base, dir0, seg_count, length_layer, scale, z)

    # convert polyline to a Rhino curve so later functions (offset etc.) work well
    crv = pl.ToNurbsCurve()

    # store main curve
    C.append(crv)
    Cmeta.append({"kind":"main","layer":layer})#new#tag as main

    # DEBUG LINES: store each polyline segment as a Line for visualization
    for i in range(pl.Count - 1):
        L.append(rg.Line(pl[i], pl[i + 1]))

    # OPTIONAL BRANCH SYSTEM
    if branch_prob > 0:
        # iterate over vertices of the polyline (excluding first and last)
        for i in range(1, pl.Count - 1):

            # branch condition
            if i % branch_interval == 0 and random.random() < branch_prob:

                # direction of the current local segment
                bdir = rg.Vector3d(pl[i + 1] - pl[i])
                bdir.Unitize()

                # force it to one of the axis directions (so branches remain 90° locked)
                if abs(bdir.X) > abs(bdir.Y):
                    bdir = rg.Vector3d(1 if bdir.X > 0 else -1, 0, 0)
                else:
                    bdir = rg.Vector3d(0, 1 if bdir.Y > 0 else -1, 0)

                # create a smaller spiral starting from that branch point
                bpl = create_spiral_polyline(
                    pl[i],
                    bdir,
                    int(seg_count * 0.25),
                    length_layer * 0.55,
                    scale,
                    z
                )

                # store branch curve + metadata
                bcrv = bpl.ToNurbsCurve()
                C.append(bcrv)
                Cmeta.append({"kind":"branch","layer":layer,"at_i":i})#new#tag as branch

                # debug lines for branch
                for j in range(bpl.Count - 1):
                    L.append(rg.Line(bpl[j], bpl[j + 1]))

    # shrink the layer’s initial length as we go up
    length_layer = max(3.0, length_layer * 0.86)


#debug#counters so we know where it dies
fail_polyline = 0
fail_outline = 0
fail_closed = 0
fail_extrude = 0

OUTLINE_DBG = []#debug#see the outline curves
OC_DBG = []#debug#see the inset curves we try to extrude

#build terraces from curves
#B will store final Breps (extruded terrace plates)
B = []

for idx, crv in enumerate(C):
    meta = Cmeta[idx]#new#metadata aligned with this curve
    is_branch = (meta.get("kind","main") == "branch")#new#now ALWAYS defined

    #Step 1) Convert curve to polyline for stable offsets
    ok, pl = crv.TryGetPolyline()#fix#TryGetPolyline has NO args in python3
    if not ok:
        fail_polyline += 1
        #fallback: manually sample the curve into points
        t0, t1 = crv.Domain.T0, crv.Domain.T1
        pts = [crv.PointAt(t0)]
        stepsN = 120
        for k in range(1, stepsN + 1):
            t = t0 + (t1 - t0) * (k / float(stepsN))
            pts.append(crv.PointAt(t))
        pl = rg.Polyline(pts)

    #Step 2) create a closed outline ribbon from polyline
    outline = offset_closed_outline(pl, terrace_w)

    #outline#guard instead of continue (GH-safe)
    if outline is None:
        fail_outline += 1
        outline_ok = False
    else:
        outline_ok = True
        OUTLINE_DBG.append(outline)#debug#visualize outlines
        if not outline.IsClosed:
            fail_closed += 1

    if outline_ok:

        #new#branches should "die out": fewer steps + faster shrink + random stop
        steps_here = terrace_steps
        shrink_here = shrink
        if is_branch:
            steps_here = max(1, int(round(terrace_steps * branch_steps_mul)))#new#branch gets fewer steps
            shrink_here = min(0.999, max(0.50, shrink * branch_shrink_mul))#new#branch shrinks faster but stays < 1

        #Step 3) make multiple inset outlines (terrace steps)
        outlines = []
        for s in range(steps_here):

            #new#random branch kill AFTER some steps (gives partial terraces that stop)
            if is_branch and s >= branch_kill_after and random.random() < branch_kill_p:
                break

            c = outline if s == 0 else scale_curve_2d(outline, (shrink_here ** s))
            if not c:
                break

            #new#area stop so tiny outlines don't become needle pillars
            amp = rg.AreaMassProperties.Compute(c)
            if (not amp) or (amp.Area < min_step_area):
                break

            outlines.append(c)
            OC_DBG.append(c)#debug#inspect inset curves

        #Step 4) extrude each outline and stack them upwards
        for s, oc in enumerate(outlines):
            if oc is None or not oc.IsClosed:
                continue

            ocT = oc.DuplicateCurve()
            ocT.Transform(rg.Transform.Translation(0, 0, s * terrace_h))

            #extrude the closed curve into a solid "plate"
            ext = rg.Extrusion.Create(ocT, terrace_h, True)

            #debug#extrusion can still fail if curve isn't valid/closed/planar
            if not ext:
                fail_extrude += 1
            else:
                B.append(ext.ToBrep())

#N = how many breps we ended with
N = len(B)
DBG = (fail_polyline, fail_outline, fail_closed, fail_extrude)
