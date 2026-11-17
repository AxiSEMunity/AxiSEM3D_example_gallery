"""Example 6 rupture: make every sub-source locatable by AxiSEM3D v2.1.0 (16 Nov 2026).

SE_Model::locateInplane looks for a point's host element only among the 30 elements with the nearest centroids (on each
MPI rank).  Next to the fine near-surface layers of the 2-s mesh, a sub-source inside a large element can have more than
30 other centroids nearer than its host's -> "Error locating source in mesh" (first seen: point2387).
Sufficient, decomposition-independent condition: the host's centroid is among the K nearest of ALL centroids, K < 30.
For every sub-source this script finds the host element (point in quad) and the rank of the host's centroid; sources
with rank >= KSAFE are moved along the straight line toward the host's centroid by the smallest fraction that brings
the rank below KSAFE.  Writes the new source file and a report.
usage: fix_source_locations.py IN.yaml OUT.yaml REPORT.txt
"""
import re, sys
import h5py, numpy as np
from scipy.spatial import cKDTree

MESH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "input", "AxiSEMCartesian_sfba_2018_2s.e")
RSURF, KSAFE, KQ = 6371100.0, 25, 400
fin, fout, frep = sys.argv[1:4]

with h5py.File(MESH, "r") as h:
    names = list(h.keys())
    if "coord" in h: xy = h["coord"][:]; X, Y = xy[0], xy[1]
    else: X, Y = h["coordx"][:], h["coordy"][:]
    con = np.concatenate([h[k][:] for k in sorted(names) if k.startswith("connect")]) - 1
print("mesh: nodes", len(X), "elements", len(con), "s", X.min(), X.max(), "z", Y.min(), Y.max(), flush=True)
P = np.stack([X[con], Y[con]], axis=-1)                       # (ne, 4, 2)
C = P.mean(axis=1); tree = cKDTree(C)


def inside(p, q, tol=1e-6):
    """point p in convex quad q (4 x 2, either orientation), boundary included"""
    d = np.roll(q, -1, axis=0) - q; r = p - q
    cr = d[:, 0] * r[:, 1] - d[:, 1] * r[:, 0]
    scale = np.abs(d).max() ** 2 * tol
    return bool(np.all(cr >= -scale) or np.all(cr <= scale))


def host_and_rank(p):
    dist, idx = tree.query(p, k=KQ)
    for rank, ie in enumerate(idx):
        if inside(p, P[ie]): return int(ie), rank
    return -1, -1


txt = open(fin).read().split("\n")
src = []                                                         # (line of distance_azimuth, line of depth, s, az, depth, name)
name = None
for i, line in enumerate(txt):
    m = re.match(r"^    - (\S+):\s*$", line)
    if m: name = m.group(1)
    m = re.match(r"^(\s*)distance_azimuth: \[([-0-9.eE+]+),\s*([-0-9.eE+]+)\]", line)
    if m:
        md = re.match(r"^(\s*)depth: ([-0-9.eE+]+)", txt[i + 1])
        src.append([i, i + 1, float(m.group(2)), float(m.group(3)), float(md.group(2)), name])
print("sources", len(src), flush=True)
rep, nfix, nout, ranks, shifts = [], 0, 0, [], []
for k, (il, idl, s, az, dep, nm) in enumerate(src):
    p = np.array([s, RSURF - dep])
    ie, rk = host_and_rank(p)
    if ie < 0:
        nout += 1; rep.append(f"{nm}: NO HOST among the {KQ} nearest (s {s:.2f} m, depth {dep:.2f} m)"); continue
    ranks.append(rk)
    if rk < KSAFE: continue
    for t in np.arange(0.02, 1.0001, 0.02):
        q = p + t * (C[ie] - p)
        ie2, rk2 = host_and_rank(q)
        if ie2 == ie and rk2 < KSAFE: break
    else:
        rep.append(f"{nm}: could not be fixed (rank {rk})"); nout += 1; continue
    sh = float(np.hypot(*(q - p))); shifts.append(sh); nfix += 1
    ind = re.match(r"^(\s*)", txt[il]).group(1)
    txt[il] = f"{ind}distance_azimuth: [{q[0]:.4f}, {az:.4f}]"
    txt[idl] = f"{ind}depth: {RSURF - q[1]:.4f}"
    rep.append(f"{nm}: host centroid rank {rk} -> {rk2}; s {s:.2f} -> {q[0]:.2f} m, depth {dep:.2f} -> {RSURF - q[1]:.2f} m; moved {sh:.1f} m "
               f"(host element {np.ptp(P[ie][:, 0]):.0f} m x {np.ptp(P[ie][:, 1]):.0f} m)")
ranks = np.array(ranks); shifts = np.array(shifts) if shifts else np.zeros(1)
head = [f"sub-sources {len(src)}; host found for {len(ranks)}; host-centroid rank: median {np.median(ranks):.0f}, max {ranks.max()}",
        f"rank >= {KSAFE} (moved): {nfix}; not fixable or no host: {nout}",
        f"shift of the moved sources: median {np.median(shifts):.1f} m, max {shifts.max():.1f} m", ""]
open(frep, "w").write("\n".join(head + rep) + "\n"); open(fout, "w").write("\n".join(txt))
print("\n".join(head)); print("\n".join(rep[:12]))
