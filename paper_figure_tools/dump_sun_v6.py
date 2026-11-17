"""Fig. 8, panel (c): snapshots of the Sun run (source half duration 382 s = mesh period).
Reads the element output (U3 at phi = 0 every 50 s, element centers) and writes the
field at the snapshot times of the figure (397, 2046, 3695, 5344, 6993, 8642 s after the origin), unfiltered: the source holds no
content below the mesh period. Also a quick-look of more times.
usage: dump_sun_v6.py XS_DIR OUT.npz QUICKLOOK.png"""
import glob, os, sys
import numpy as np, h5py
XS, OUT, QL = sys.argv[1:4]
GC = 4
T_SNAP = np.array([397.2, 2046.1, 3694.9, 5343.8, 6992.7, 8641.5])
fs = sorted(glob.glob(os.path.join(XS, 'axisem3d_synthetics.nc.rank*')))
S, Z, Wt, t = [], [], [], None
for f in fs:
    with h5py.File(f, 'r') as d:
        if t is None: t = d['data_time'][:].astype(float)
        co = d['list_element_coords'][:]; w = d['data_wave__NaG=1'][:, 0, GC, 0, :]
        S.append(co[:, GC, 0]); Z.append(co[:, GC, 1]); Wt.append(w.astype(np.float32))
S = np.concatenate(S); Z = np.concatenate(Z); Wt = np.concatenate(Wt, axis=0)
idx = np.array([int(np.argmin(np.abs(t - x))) for x in T_SNAP])
print(f'{len(fs)} files, {S.size} element centers, {t.size} samples, t = {t[0]:.1f}..{t[-1]:.1f} s; snapshot times (s):', np.round(t[idx], 1))
assert np.all(np.abs(t[idx] - T_SNAP) < 30), 'snapshot times not on the 50-s grid'
assert np.isfinite(Wt).all(), 'NaN in the output'
np.savez_compressed(OUT, S=S, Z=Z, W=Wt[:, idx], t_snap=t[idx], t_all=t)
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
qi = [int(np.argmin(np.abs(t - x))) for x in np.arange(400, t[-1], 1200)][:10]
fig, axs = plt.subplots(2, 5, figsize=(17, 7)); R = 6.96e8
for ax, i in zip(axs.ravel(), qi):
    v = Wt[:, i].astype(np.float64); vm = max(float(np.percentile(np.abs(v), 98)), 1e-30)
    ax.scatter(np.r_[S, -S] / 1e8, np.r_[Z, Z] / 1e8, c=np.r_[v, v], s=0.3, cmap='RdBu_r', vmin=-vm, vmax=vm, rasterized=True, linewidths=0, marker='s')
    th = np.linspace(0, 2 * np.pi, 300); ax.plot(R / 1e8 * np.cos(th), R / 1e8 * np.sin(th), 'k', lw=0.5)
    ax.set_aspect('equal'); ax.axis('off'); ax.set_title(f't = {t[i] / 3600:.2f} h', fontsize=9)
fig.savefig(QL, dpi=110, bbox_inches='tight'); print('wrote', OUT, QL)
