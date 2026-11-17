"""Example 6: build the rupture source list for the 2-s run from the shipped 36,246-sub-source list (16 Nov 2026).

Two things are done, each with a check:

1. FRAME.  AxiSEM3D reads the moment tensor of a source given by `distance_azimuth` as [Mzz, Mss, Mpp, Mzs, Mzp, Msp] in
   the LOCAL cylindrical frame at the source's azimuth (Source::computeQzsp applies no rotation for source-centered
   locations; Mechanism::release places the source with exp(-i alpha phi)).  The shipped list holds tensors in ONE FIXED
   geographic frame (checked: the horizontal principal axes are the same for all sources, also for
   those near the axis whose azimuths cover the full circle).  Read as (vertical, north, east) they are right-lateral
   slip on a plane striking about N38W, i.e. on the fault traced by the source positions (x = east, y = north,
   azimuth counter-clockwise from east; frame verified with summits).  Here every tensor is rotated into
   the local frame the solver expects.  The sense of the vertical axis (up or down) is chosen by comparing the nodal
   plane of the summed tensor with the plane fitted to the source positions.

2. NUMBER OF SUB-SOURCES.  The solver applies every source at every time step at the cost of one element evaluation;
   36,246 sub-sources made the 2-s run 16 times slower (67 h).  The sub-sources (about 110 m apart, designed for 1-s runs)
   are merged on a grid of --cell metres along strike and depth: tensors are summed in the fixed frame (total moment
   tensor preserved exactly), position and rupture time are moment-weighted means.

usage: build_rupture_2s.py IN.yaml OUT.yaml OUT_TABLE.csv [--cell 400]
"""
import argparse, re
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("fin"); ap.add_argument("fout"); ap.add_argument("ftab")
ap.add_argument("--cell", type=float, default=400.0)
a = ap.parse_args()
txt = open(a.fin).read()
head = txt[:txt.index("list_of_sources:")]
da = np.array(re.findall(r"distance_azimuth: \[([-0-9.eE+]+),\s*([-0-9.eE+]+)\]", txt), float)
dep = np.array(re.findall(r"distance_azimuth: \[[^\]]*\]\s*\n\s*depth: ([-0-9.eE+]+)", txt), float)
M = np.array([[float(v) for v in m.split(",")] for m in re.findall(r"data: \[([^\]]+)\]", txt)])
tsh = np.array(re.findall(r"time_shift: ([-0-9.eE+]+)", txt), float)
unit = set(re.findall(r"\n\s*unit: ([-0-9.eE+]+)", txt)); hd = set(re.findall(r"half_duration: ([-0-9.eE+]+)", txt))
assert len(da) == len(dep) == len(M) == len(tsh) and len(unit) == 1 and len(hd) == 1, (len(da), len(dep), len(M), len(tsh), unit, hd)
unit = float(unit.pop()); n0 = len(da)
s, phi = da[:, 0], da[:, 1]
E, N = s * np.cos(phi), s * np.sin(phi)                                   # east, north (m)


def full(m, zsign):
    """6 components [11, 22, 33, 12, 13, 23] with 1 = vertical, 2 = north, 3 = east -> 3x3 in (east, north, up)"""
    T = np.zeros(m.shape[:-1] + (3, 3))
    T[..., 2, 2] = m[..., 0]; T[..., 1, 1] = m[..., 1]; T[..., 0, 0] = m[..., 2]
    T[..., 2, 1] = T[..., 1, 2] = zsign * m[..., 3]; T[..., 2, 0] = T[..., 0, 2] = zsign * m[..., 4]
    T[..., 1, 0] = T[..., 0, 1] = m[..., 5]
    return T


def planes(T):
    """the two nodal planes of the best double couple: (normal, slip) pairs in (east, north, up)"""
    w, v = np.linalg.eigh(T - np.trace(T) / 3 * np.eye(3)); t, p = v[:, 2], v[:, 0]
    return [((t + p) / np.sqrt(2), (t - p) / np.sqrt(2)), ((t - p) / np.sqrt(2), (t + p) / np.sqrt(2))]


def sdr(n, d):
    """strike, dip, rake (deg, Aki & Richards) from fault normal n and slip d in (east, north, up)"""
    if n[2] > 0: n, d = -n, -d                                           # normal pointing down = hanging wall side
    dip = np.degrees(np.arccos(-n[2])); strike = np.degrees(np.arctan2(-n[1], n[0])) % 360
    sv = np.array([np.sin(np.radians(strike)), np.cos(np.radians(strike)), 0.0])
    dv = np.cross(-n, sv); rake = np.degrees(np.arctan2(np.dot(d, dv), np.dot(d, sv)))
    return strike, dip, rake


# ---- plane through the source positions
X = np.c_[E, N, -dep]; c0 = X.mean(axis=0); u_, sv_, vt = np.linalg.svd(X - c0, full_matrices=False)
ngeo = vt[2] / np.linalg.norm(vt[2])
print(f"sub-sources {n0}; plane through the positions: strike {sdr(ngeo, vt[0])[0]:.1f} or {(sdr(ngeo, vt[0])[0]+180)%360:.1f} deg, "
      f"dip {sdr(ngeo, vt[0])[1]:.1f} deg; rms distance from the plane {sv_[2]/np.sqrt(n0):.0f} m")
best = None
for zsign, lab in ((1.0, "vertical axis up"), (-1.0, "vertical axis down")):
    Tsum = full(M, zsign).sum(axis=0)
    for n_, d_ in planes(Tsum):
        agree = abs(float(np.dot(n_, ngeo))); st, dp, rk = sdr(n_, d_)
        print(f"  {lab:20s}: nodal plane strike {st:6.1f}, dip {dp:5.1f}, rake {rk:7.1f} | |cos| to the geometric plane {agree:.4f}")
        if best is None or agree > best[0]: best = (agree, zsign, lab, st, dp, rk)
agree, zsign, lab, st, dp, rk = best
print(f"chosen: {lab} (nodal plane within {np.degrees(np.arccos(min(agree, 1.0))):.1f} deg of the geometric plane; strike {st:.1f}, dip {dp:.1f}, rake {rk:.1f})")
T = full(M, zsign)                                                        # (n0, 3, 3) in (east, north, up)
m0 = np.sqrt(0.5 * (T ** 2).sum(axis=(1, 2)))

# ---- merge on a grid along strike and depth
ev = vt[0] if abs(vt[0][2]) < abs(vt[1][2]) else vt[1]                   # along-strike unit vector (the more horizontal one)
ev = np.array([ev[0], ev[1], 0.0]); ev /= np.linalg.norm(ev)
u = (X - c0) @ ev
iu = np.floor((u - u.min()) / a.cell).astype(int); idp = np.floor(dep / a.cell).astype(int)
key = iu * 100000 + idp; uk, inv = np.unique(key, return_inverse=True); nm = len(uk)
W = np.bincount(inv, weights=m0, minlength=nm)
Tm = np.zeros((nm, 3, 3)); np.add.at(Tm, inv, T)
pos = np.stack([np.bincount(inv, weights=m0 * q, minlength=nm) / W for q in (E, N, dep)], axis=1)
tm = np.bincount(inv, weights=m0 * tsh, minlength=nm) / W
tspread = np.sqrt(np.maximum(np.bincount(inv, weights=m0 * tsh ** 2, minlength=nm) / W - tm ** 2, 0.0))
cnt = np.bincount(inv, minlength=nm)
sm, pm = np.hypot(pos[:, 0], pos[:, 1]), np.arctan2(pos[:, 1], pos[:, 0]) % (2 * np.pi)
# ---- rotate into the local frame (z up, s, phi) at the merged source
es = np.stack([np.cos(pm), np.sin(pm), np.zeros(nm)], axis=1); ep = np.stack([-np.sin(pm), np.cos(pm), np.zeros(nm)], axis=1)
ez = np.tile([0.0, 0.0, 1.0], (nm, 1))
comp = lambda p_, q_: np.einsum("ni,nij,nj->n", p_, Tm, q_)
L = np.stack([comp(ez, ez), comp(es, es), comp(ep, ep), comp(ez, es), comp(ez, ep), comp(es, ep)], axis=1)

mw = lambda x_: (np.log10(x_ * unit) - 9.1) / 1.5
sc = lambda Tx: np.sqrt(0.5 * (Tx ** 2).sum())
print(f"merged on a {a.cell:.0f}-m grid: {n0} -> {nm} sub-sources (1 to {cnt.max()} per cell, median {np.median(cnt):.0f})")
print(f"summed tensor: before {sc(T.sum(axis=0)):.6e}, after {sc(Tm.sum(axis=0)):.6e} (file units) -> M0 = {sc(Tm.sum(axis=0))*unit:.3e} N m, Mw {mw(sc(Tm.sum(axis=0))):.2f}")
print(f"sum of scalar moments: before {m0.sum()*unit:.3e} N m (Mw {mw(m0.sum()):.2f}), after {np.sqrt(0.5*(Tm**2).sum(axis=(1,2))).sum()*unit:.3e} N m")
print(f"rupture time within a cell: rms spread median {np.median(tspread):.3f} s, max {tspread.max():.3f} s; rupture times {tm.min():.2f}-{tm.max():.2f} s")
far = sm > 3000
hs = np.sqrt(0.5 * (L[:, 1] ** 2 + L[:, 2] ** 2) + L[:, 5] ** 2)
print(f"local frame, sources beyond 3 km: |Msp| / horizontal norm = {np.median(np.abs(L[far, 5]) / hs[far]):.3f} (1 = pure strike-slip on a radial plane); "
      f"before the rotation the same ratio was {np.median(np.abs(M[s > 3000, 5]) / np.sqrt(0.5*(M[s > 3000, 1]**2 + M[s > 3000, 2]**2) + M[s > 3000, 5]**2)):.3f}")

# ---- write
hdv = hd.pop(); o = np.lexsort((pos[:, 2], u.min() + (uk // 100000) * a.cell))
blocks = []
for j, k in enumerate(o):
    blocks.append(f"""    - sub{j + 1:05d}:
        location:
            distance_azimuth: [{sm[k]:.4f}, {pm[k]:.6f}]
            depth: {pos[k, 2]:.4f}
            ellipticity: false
            depth_below_solid_surface: true
            undulated_geometry: false
        mechanism:
            type: MOMENT_TENSOR
            data: [{', '.join(f'{v:.6e}' for v in L[k])}]
            unit: {unit:g}
        source_time_function:
            class_name: GaussianSTF
            half_duration: {hdv}
            decay_factor: 1.628
            time_shift: {tm[k]:.4f}
            use_derivative_integral: ERF""")
note = ("# Source list of the Example 6 rupture run (Fig. 6), built by build_rupture_2s.py from the\n"
        f"# {n0}-sub-source list shipped with AxiSEM3D v2.1.0: tensors rotated from the fixed (vertical, north, east) frame of\n"
        f"# that list into the local (z, s, phi) frame the solver reads, then merged on a {a.cell:.0f}-m grid ({nm} sub-sources).\n"
        f"# Moments in dyne cm (unit {unit:g}).  Total moment {sc(Tm.sum(axis=0))*unit:.3e} N m, Mw {mw(sc(Tm.sum(axis=0))):.2f}.\n")
open(a.fout, "w").write(note + head + "list_of_sources:\n" + "\n".join(blocks) + "\n")
np.savetxt(a.ftab, np.c_[pos[o, 0], pos[o, 1], pos[o, 2], sm[o], pm[o], tm[o], np.sqrt(0.5 * (Tm[o] ** 2).sum(axis=(1, 2))) * unit, cnt[o]],
           header="east_m north_m depth_m distance_m azimuth_rad rupture_time_s moment_Nm n_merged", fmt="%.4f %.4f %.4f %.4f %.6f %.4f %.6e %d")
print("wrote", a.fout, "and", a.ftab)
