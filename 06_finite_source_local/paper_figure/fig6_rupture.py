"""Fig. 6 (Example 6): ground velocity of the Hayward-fault rupture in the San Francisco Bay Area model.
Reads the station output of the rupture run (displacement, ENZ, 250-m receiver grid on a 45-km disk), differentiates to
velocity and draws the manuscript's four panels: (a) vertical, (b) east, (c) north velocity at one instant, (d) peak
ground velocity (maximum of the three-component velocity magnitude over the record).
Frame: x = east, y = north of the model center (37.85517, -122.25684) (verified with summits).
Shoreline: GSHHG full resolution, level 1 (no lakes, no rivers), azimuthal equidistant about the same center.
No transparency anywhere; PDF with TrueType text.
usage: HDF5_USE_FILE_LOCKING=FALSE python fig6_rupture.py RUN_DIR OUT_STEM [--tsnap 16] [--pgvtmax 16] [--vmax V] [--pgvmax P] [--partial]
       (Fig. 6: --tsnap 16 --pgvtmax 16)
"""
import argparse, glob, os, re
import h5py, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

ap = argparse.ArgumentParser(); ap.add_argument("run"); ap.add_argument("out")
ap.add_argument("--tsnap", type=float, default=16.0); ap.add_argument("--vmax", type=float, default=None)
ap.add_argument("--pgvmax", type=float, default=None); ap.add_argument("--partial", action="store_true")
ap.add_argument("--pgvtmax", type=float, default=None, help="peak ground velocity over t <= PGVTMAX only")
ap.add_argument("--group", default="sfba_local"); ap.add_argument("--dx", type=float, default=250.0)
a = ap.parse_args()
HERE = os.path.dirname(os.path.abspath(__file__))
LAT0, LON0, RE, DX, RMAX = 37.85517, -122.25684, 6371.0e3, a.dx, 45000.0
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 8,
                     "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.linewidth": 0.6,
                     "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.transparent": False})


def aeqd(lon, lat):
    la, lo, la0, lo0 = map(np.radians, (lat, lon, LAT0, LON0))
    c = np.arccos(np.clip(np.sin(la0) * np.sin(la) + np.cos(la0) * np.cos(la) * np.cos(lo - lo0), -1, 1))
    az = np.arctan2(np.sin(lo - lo0) * np.cos(la), np.cos(la0) * np.sin(la) - np.sin(la0) * np.cos(la) * np.cos(lo - lo0))
    return RE * c * np.sin(az) / 1e3, RE * c * np.cos(az) / 1e3


# ---- stations and data
pos = {}
for line in open(os.path.join(a.run, "input", "STATIONS_OUTPUT.txt")):
    if line.startswith("#"): continue
    p = line.split(); pos[f"{p[1]}.{p[0]}"] = (float(p[2]), float(p[3]))
N = int(round(2 * RMAX / DX)) + 1
names, chunks, t, chan = [], [], None, None
for fn in sorted(glob.glob(os.path.join(a.run, "output", "stations", a.group, "axisem3d_synthetics.nc.rank*"))):
    try:
        with h5py.File(fn, "r") as f:
            st = [b"".join(s).decode("utf-8", "ignore").strip().strip("\x00") for s in f["list_station"][:]]
            if not st: continue
            tt = f["data_time"][:]
            if a.partial: tt = tt[tt < 1e30]
            d = f["data_wave"][:, :, :len(tt)]
            if chan is None: chan = [b"".join(c).decode().strip().strip("\x00") for c in f["list_channel"][:]]
    except OSError as e:
        print("skip", os.path.basename(fn), e); continue
    names += st; chunks.append(d.astype(np.float32)); t = tt if t is None or len(tt) < len(t) else t
nt = min(c.shape[2] for c in chunks); t = np.asarray(t[:nt], float)
U = np.concatenate([c[:, :, :nt] for c in chunks], axis=0)
print("stations read", len(names), "of", len(pos), "| channels", chan, "| samples", nt, "| t", t[0], t[-1], flush=True)
iE, iN, iZ = [next(i for i, c in enumerate(chan) if c.upper().endswith(k)) for k in ("1", "2", "3")] if not any("E" in c.upper() for c in chan) \
    else [next(i for i, c in enumerate(chan) if c.upper().endswith(k)) for k in ("E", "N", "Z")]
dt = t[1] - t[0]
V = np.gradient(U, dt, axis=2)                                   # velocity (m/s)
speed = np.sqrt((V ** 2).sum(axis=1))
pgv = speed[:, t <= a.pgvtmax].max(axis=1) if a.pgvtmax else speed.max(axis=1)
ks = int(np.argmin(np.abs(t - a.tsnap)))
xy = np.array([pos[n] for n in names]); ix = np.rint((xy[:, 0] + RMAX) / DX).astype(int); iy = np.rint((xy[:, 1] + RMAX) / DX).astype(int)


def grid(v):
    g = np.full((N, N), np.nan, dtype=np.float32); g[iy, ix] = v; return g


print(f"t_snap {t[ks]:.2f} s | max |v| at snapshot {np.abs(V[:, :, ks]).max():.3g} m/s | PGV max {pgv.max():.3g}, 99th pct {np.percentile(pgv, 99):.3g}, "
      f"median {np.median(pgv):.3g} m/s", flush=True)

# ---- fault trace and hypocenter from the source list (x = s cos(az), y = s sin(az))
src = open(os.path.join(a.run, "input", "inparam.source.yaml")).read()
da = np.array(re.findall(r"distance_azimuth: \[([-0-9.eE+]+),\s*([-0-9.eE+]+)\]", src), float)
HAVE_FAULT = len(da) > 1
if not HAVE_FAULT: da = np.zeros((2, 2))
dep = np.array(re.findall(r"distance_azimuth: \[[^\]]*\]\s*\n\s*depth: ([-0-9.eE+]+)", src), float) if HAVE_FAULT else np.zeros(2)
tsh = np.array(re.findall(r"time_shift: ([-0-9.eE+]+)", src), float)[:len(da)] if HAVE_FAULT else np.zeros(2)
fx, fy = da[:, 0] * np.cos(da[:, 1]) / 1e3, da[:, 0] * np.sin(da[:, 1]) / 1e3
top = (dep < np.percentile(dep, 3)) if HAVE_FAULT else np.ones(len(dep), bool)
along = fx[top] * np.cos(np.radians(-55.0)) + fy[top] * np.sin(np.radians(-55.0)); o = np.argsort(along)
hyp = int(np.argmin(tsh))
if HAVE_FAULT: print(f"fault: {len(da)} sub-sources, surface trace from ({fx[top][o][0]:.1f}, {fy[top][o][0]:.1f}) to ({fx[top][o][-1]:.1f}, {fy[top][o][-1]:.1f}) km; "
      f"hypocenter at ({fx[hyp]:.1f}, {fy[hyp]:.1f}) km, depth {dep[hyp]/1e3:.1f} km", flush=True)

# ---- shoreline
seg, cur = [], []
for line in open(os.path.join(HERE, "coast_sfba_gshhg_full_L1.txt")):
    if line.startswith(">"):
        if cur: seg.append(np.array(cur)); cur = []
    else:
        cur.append([float(v) for v in line.split()[:2]])
if cur: seg.append(np.array(cur))

# ---- figure
vmax = a.vmax or float(np.percentile(np.abs(V[:, :, ks]), 99.7))
pmax = a.pgvmax or float(np.percentile(pgv, 99.5))
seq = LinearSegmentedColormap.from_list("pgv", ["#fffdf0", "#fde9a2", "#fbb74f", "#ee6a2c", "#c1272d", "#6a0d1f"])
div = plt.get_cmap("RdBu_r").copy(); div.set_bad("white"); seq.set_bad("white")
L = 47.0
fig, axs = plt.subplots(2, 2, figsize=(6.85, 7.2))
fig.subplots_adjust(left=0.09, right=0.99, top=0.995, bottom=0.115, wspace=0.03, hspace=0.03)
panels = (("(a) Vertical", grid(V[:, iZ, ks]), div, -vmax, vmax), ("(b) East", grid(V[:, iE, ks]), div, -vmax, vmax),
          ("(c) North", grid(V[:, iN, ks]), div, -vmax, vmax), ("(d) Peak ground velocity", grid(pgv), seq, 0.0, pmax))
lons, lats = np.array([-122.6, -122.2, -121.8]), np.array([37.5, 37.7, 37.9, 38.1])
ims = []
for ax, (lab, g, cm, lo, hi) in zip(axs.ravel(), panels):
    ims.append(ax.imshow(g, extent=(-RMAX / 1e3 - DX / 2e3, RMAX / 1e3 + DX / 2e3, -RMAX / 1e3 - DX / 2e3, RMAX / 1e3 + DX / 2e3), origin="lower",
                         cmap=cm, vmin=lo, vmax=hi, interpolation="bilinear", zorder=1))
    for lo_ in lons:
        e, n = aeqd(np.full(50, lo_), np.linspace(37.3, 38.5, 50)); ax.plot(e, n, color="0.72", lw=0.35, zorder=2)
    for la_ in lats:
        e, n = aeqd(np.linspace(-123.0, -121.4, 50), np.full(50, la_)); ax.plot(e, n, color="0.72", lw=0.35, zorder=2)
    for s_ in seg:
        e, n = aeqd(s_[:, 0], s_[:, 1]); ax.plot(e, n, color="0.15", lw=0.7, zorder=3)
    if HAVE_FAULT:
        ax.plot(fx[top][o], fy[top][o], color="k", lw=0.9, zorder=4)
        ax.plot([fx[hyp]], [fy[hyp]], marker="*", ms=7, mfc="white", mec="k", mew=0.7, zorder=5)
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.text(0.02, 0.975, lab, transform=ax.transAxes, va="top", ha="left", fontsize=8.5, zorder=6,
            bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
    xt = [aeqd(lo_, LAT0 - L / 111.19)[0] for lo_ in lons]; yt = [aeqd(LON0 - L / 87.8, la_)[1] for la_ in lats]
    ax.set_xticks(xt); ax.set_yticks(yt)
    ax.set_xticklabels([f"{-lo_:.1f}°W" for lo_ in lons]); ax.set_yticklabels([f"{la_:.2f}°N" for la_ in lats])
    ax.tick_params(length=2.5, width=0.5)
for ax in axs[0]: ax.set_xticklabels([])
for ax in axs[:, 1]: ax.set_yticklabels([])
cax1 = fig.add_axes([0.14, 0.066, 0.34, 0.014]); cax2 = fig.add_axes([0.59, 0.066, 0.34, 0.014])   # just below the tick labels
cb1 = fig.colorbar(ims[0], cax=cax1, orientation="horizontal", extend="both"); cb1.set_label("Ground velocity (m/s)", labelpad=2)
cb2 = fig.colorbar(ims[3], cax=cax2, orientation="horizontal", extend="max"); cb2.set_label("Peak ground velocity (m/s)", labelpad=2)
from matplotlib.ticker import MaxNLocator
for cb in (cb1, cb2):
    cb.locator = MaxNLocator(nbins=5); cb.update_ticks()
    cb.outline.set_linewidth(0.5); cb.ax.tick_params(length=2.5, width=0.5, labelsize=7)
fig.savefig(a.out + ".png", dpi=600, bbox_inches="tight", pad_inches=0.02); fig.savefig(a.out + ".pdf", bbox_inches="tight", pad_inches=0.02)
from PIL import Image
_im = Image.open(a.out + ".png"); _a = _im.getchannel("A") if _im.mode == "RGBA" else None
assert _a is None or _a.getextrema() == (255, 255), "non-opaque pixels"
_im = _im.convert("RGB"); _im.save(a.out + ".png", dpi=(600, 600)); _im.save(a.out + ".jpg", quality=95, subsampling=0, dpi=(600, 600))
print("wrote", a.out + ".png", "| snapshot at", round(float(t[ks]), 2), "s | vmax", vmax, "| pgv max", pmax)
