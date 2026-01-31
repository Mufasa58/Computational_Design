import Rhino.Geometry as rg
import random
import math

class OrganicGrowthSystem:
    def __init__(self, start_curve, params):
        self.current_curve = start_curve.Duplicate()
        self.iterations = int(params.get('iterations', 10))
        self.z_dist = params.get('z_dist', 1.0)
        self.prob = params.get('prob', 0.5)
        self.shrink = params.get('shrink', 0.8)
        self.size_proportion = params.get('size_proportion', 1.5)
        
        # FIXED LENGTH LOGIC
        self.fixed_L1 = params.get('fixed_L1', 2.0)
        self.fixed_L3 = params.get('fixed_L3', 1.0)
        self.jitter = params.get('jitter', 0.1) 
        
        # GRADUAL & CONSTRAINTS
        self.base_offset = params.get('xy_off', 1.0)
        self.taper_factor = params.get('taper_factor', 0.95)
        self.decay_rate = params.get('decay_rate', 0.02)
        
        # MINIMUM SIZE LOGIC (Safety Stops)
        self.initial_length = self.current_curve.GetLength()
        self.min_length = self.initial_length * 0.35 
        self.stop_threshold = 0.01 # Stop when total length is very close to zero
        
        self.current_trim = params.get('step_base', 0.1)
        self.history = []
        random.seed(params.get('seed', 1))

    def _get_random_fixed(self, base_val, iter_idx):
        decayed = base_val * (self.taper_factor ** iter_idx)
        variation = decayed * self.jitter
        return decayed + random.uniform(-variation, variation)

    def _grow_three_relative(self, start_pt, start_dir, anterior_len, is_clockwise, iter_idx):
        lines = []
        curr_pos = start_pt
        curr_dir = rg.Vector3d(start_dir)
        if not curr_dir.Unitize(): return lines

        l1 = self._get_random_fixed(self.fixed_L1, iter_idx)
        l2 = anterior_len * self.shrink
        l3 = self._get_random_fixed(self.fixed_L3, iter_idx)
        
        lengths = [max(0.001, l1), max(0.001, l2), max(0.001, l3)]
        angle = -math.pi / 2 if is_clockwise else math.pi / 2

        for i in range(3):
            curr_dir.Rotate(angle, rg.Vector3d.ZAxis)
            seg_len = lengths[i]
            new_p = curr_pos + (curr_dir * seg_len)
            lines.append(rg.LineCurve(curr_pos, new_p))
            curr_pos = new_p
            
        return lines

    def run(self):
        for i in range(self.iterations):
            self.current_curve.Domain = rg.Interval(0, 1)
            total_len = self.current_curve.GetLength()

            # 1. SAFETY STOP: Exit if curve is near zero
            if total_len < self.stop_threshold:
                break

            potential_trim = self.current_trim + (total_len * self.decay_rate)
            resulting_len = total_len - (2 * potential_trim)

            if resulting_len > self.min_length:
                self.current_trim = potential_trim
            else:
                self.current_trim = (total_len - self.min_length) / 2
            
            self.current_trim = max(0, self.current_trim)

            # Avoid trimming the curve into non-existence
            if (self.current_trim * 2) >= total_len:
                break

            t_start = self.current_trim / total_len
            t_end = 1.0 - (self.current_trim / total_len)

            shrunk_curve = self.current_curve.Trim(t_start, t_end)
            if not shrunk_curve: break

            # 2. BRANCHING
            to_join = [shrunk_curve]
            success, polyline = shrunk_curve.TryGetPolyline()
            
            if success and polyline.Count >= 2:
                p0, p_last = polyline[0], polyline[polyline.Count - 1]
                
                branch_len_start = p0.DistanceTo(polyline[1]) / self.size_proportion
                branch_len_end = p_last.DistanceTo(polyline[polyline.Count-2]) / self.size_proportion

                if branch_len_start > self.stop_threshold and random.random() < self.prob:
                    branches = self._grow_three_relative(p0, p0 - polyline[1], branch_len_start, True, i)
                    for seg in branches: seg.Reverse()
                    branches.reverse()
                    to_join = branches + to_join

                if branch_len_end > self.stop_threshold and random.random() < self.prob:
                    branches = self._grow_three_relative(p_last, p_last - polyline[polyline.Count-2], branch_len_end, False, i)
                    to_join.extend(branches)

            # 3. JOIN AND DYNAMIC OFFSET
            # FIX: Used rg.Curve.JoinCurves instead of rg.JoinCurves
            joined = rg.Curve.JoinCurves(to_join, 0.01)
            if not joined or len(joined) == 0: break
            
            temp_curve = joined[0]
            current_offset = self.base_offset * (self.taper_factor ** i)
            
            # Stop if offset becomes too small to solve
            if current_offset < 0.001: break

            offset_result = temp_curve.Offset(rg.Plane.WorldXY, current_offset, 0.01, rg.CurveOffsetCornerStyle.Sharp)
            
            if offset_result:
                self.current_curve = offset_result[0]
            else:
                break # Stop if offset fails

            # 4. RECORD
            self.current_curve.Transform(rg.Transform.Translation(0, 0, self.z_dist))
            self.history.append(self.current_curve.DuplicateCurve())

        return self.history

# =========================================================
# GH SCRIPT AREA
# =========================================================
params = {
    'iterations': N,
    'step_base': Step,
    'xy_off': XY_Off,
    'z_dist': Z_Dist,
    'seed': Seed,
    'prob': prob,
    'shrink': shrink,
    'size_proportion': size_proportion,
    'decay_rate': 0.01,
    'taper_factor': 0.98,
    'fixed_L1': 1,
    'fixed_L3': 1,
    'jitter': 0.05
}

system = OrganicGrowthSystem(C, params)
a = system.run()