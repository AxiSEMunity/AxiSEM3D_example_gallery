"""Fig. 7 (Example 7): P-wave traveltime kernel for v_p, and its test against direct simulation.
(a) 1-D PREM, (b) 3-D model (S362ANI x 5 + Moho undulation): relative v_p kernel on the source-receiver plane,
    K = 2 rho vp^2 K_lambda / tau (s m^-3; dT = int K dln(vp) dV); dashed: the test perturbation (half-maximum contour);
(c) the 3-D v_p model on the same plane (percent from PREM);
(d) traveltime change of P predicted by the kernel against the change simulated directly, for the four perturbations
    (+-0.5 %, +-1 %) in both models.
Data: data/kvp_<dim>_slice_part*.npz (eval_kvp.py), data/test_kernels_<dim>.npz (test_kernels.py),
S362ANI_radius.nc of the example.   usage: plot_fig7.py OUT_STEM
"""
import glob, os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import h5py

HERE = os.path.dirname(os.path.abspath(__file__)); D = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'data')
OUT = sys.argv[1]
S362 = os.path.join(HERE, '..', 'kernel_3D', 'input_forward', 'S362ANI_radius.nc')
R_E, TAU = 6371.0, 2.0
CLAT, CLON, SH, RC, SR = 0.0, 20.0, 5.0, 6371.0 - 960.0, 250.0
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': 8,
                     'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7, 'axes.linewidth': 0.6,
                     'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.transparent': False, 'image.interpolation': 'nearest'})


def load_slice(dim):
    fns = sorted(glob.glob(f'{D}/kvp_{dim}_slice_part*.npz')); z0 = np.load(fns[0]); shape = tuple(int(x) for x in z0['shape'])
    n = int(np.prod(shape)); A = {k: np.full(n, np.nan) for k in ('K_lambda', 'rho1d', 'vp1d', 'r', 'lon')}
    for fn in fns:
        z = np.load(fn)
        for k in A: A[k][z['index']] = z[k]
    A = {k: v.reshape(shape)[:, 0, :] for k, v in A.items()}
    return A, float(z0['TAU'])


def s362_dvp(r_km, lon):
    from scipy.interpolate import RegularGridInterpolator
    try:
        with h5py.File(S362, 'r') as d:
            rr = d['radius'][()].astype(float); la = d['latitude'][()].astype(float); lo = d['longitude'][()].astype(float); dvs = d['dvs'][()].astype(float)
    except OSError:
        from scipy.io import netcdf_file
        with netcdf_file(S362, 'r', mmap=False) as d:
            rr = d.variables['radius'][:].astype(float); la = d.variables['latitude'][:].astype(float); lo = d.variables['longitude'][:].astype(float); dvs = d.variables['dvs'][:].astype(float)
    o = np.argsort(rr); f = RegularGridInterpolator((rr[o], la, lo), dvs[o], bounds_error=False, fill_value=0.0)
    return 5.0 * f(np.c_[r_km.ravel() * 1e3, np.zeros(r_km.size), lon.ravel()]).reshape(r_km.shape)


def section_xy(r_km, lon):
    th = np.radians(lon - CLON); return r_km * np.sin(th), r_km * np.cos(th)


K, dvp3 = {}, None
for dim in ('1d', '3d'):
    A, tau = load_slice(dim); r_km = A['r'] / 1e3
    vp = A['vp1d']
    if dim == '3d':
        dvp3 = s362_dvp(r_km, A['lon']); vp = vp * (1 + dvp3)
    K[dim] = (2 * A['rho1d'] * vp ** 2 * A['K_lambda'] / tau, r_km, A['lon'])
far = lambda lon, r: (np.minimum(np.abs(lon - 0.0), np.abs(lon - 40.0)) > 6.0) | (r < R_E - 400.0)

vmax = 0.8 * float(max(np.nanmax(np.abs(k[far(lon, r)])) for k, r, lon in K.values()))
print(f'kernel color limit +-{vmax:.3e} s/m^3 (0.8 x the largest |K| away from source and receiver)')
RMIN, HALF = R_E - 1900.0, 25.0
# layout in inches: each section panel is exactly as tall as its equal-aspect drawing (no slack), and each color bar
# sits 0.07 in below its section
FW, SW = 6.85, 0.495 * 6.85                                          # figure width, section panel width
SECH = SW * ((R_E + 160) - (RMIN * np.cos(np.radians(HALF)) - 40)) / (2 * R_E * np.sin(np.radians(HALF)) + 120)
CBH, GAP, LAB = 0.065, 0.07, 0.38                                    # color bar height, bar-section gap, room for bar labels
Y_CB2 = 0.36; Y_C = Y_CB2 + CBH + GAP; Y_CB1 = Y_C + SECH + LAB; Y_AB = Y_CB1 + CBH + GAP; FH = Y_AB + SECH + 0.05
inch = lambda x, y, w, h: [x / FW, y / FH, w / FW, h / FH]
fig = plt.figure(figsize=(FW, FH))
axs = [fig.add_axes(inch(0.0, Y_AB, SW, SECH)), fig.add_axes(inch(0.505 * FW, Y_AB, SW, SECH)), fig.add_axes(inch(0.0, Y_C, SW, SECH))]
for a_ in axs: a_.set_anchor('S')
th_b = np.radians(np.linspace(-HALF, HALF, 300))


def frame(ax):
    for rr_ in (R_E, RMIN):
        ax.plot(rr_ * np.sin(th_b), rr_ * np.cos(th_b), color='0.25', lw=0.7)
    for t_ in (-HALF, HALF):
        ax.plot([RMIN * np.sin(np.radians(t_)), R_E * np.sin(np.radians(t_))], [RMIN * np.cos(np.radians(t_)), R_E * np.cos(np.radians(t_))], color='0.25', lw=0.7)
    ax.set_aspect('equal'); ax.axis('off')
    ax.set_xlim(R_E * np.sin(np.radians(-HALF)) - 60, R_E * np.sin(np.radians(HALF)) + 60); ax.set_ylim(RMIN * np.cos(np.radians(HALF)) - 40, R_E + 160)


def sector(r_km, lon):
    rows = r_km[:, 0] >= RMIN - 5; cols = np.abs(lon[0] - CLON) <= HALF + 0.05
    return np.ix_(rows, cols)


for ax, key, lab in ((axs[0], '1d', '(a)'), (axs[1], '3d', '(b)')):
    k, r_km, lon = K[key]; X, Y = section_xy(r_km, lon); s_ = sector(r_km, lon)
    im = ax.pcolormesh(X[s_], Y[s_], k[s_] / 1e-16, cmap='RdBu_r', norm=TwoSlopeNorm(0, -vmax / 1e-16, vmax / 1e-16), shading='gouraud', rasterized=True)
    TT, RR = np.meshgrid(np.linspace(0, 40, 400), np.linspace(RMIN, R_E, 300))
    g = np.exp(-0.5 * ((TT - CLON) / SH) ** 2) * np.exp(-0.5 * ((RR - RC) / SR) ** 2)
    Xb, Yb = section_xy(RR, TT); ax.contour(Xb, Yb, g, levels=[0.5], colors='k', linewidths=0.7, linestyles='--')
    frame(ax)
    xs, ys = section_xy(np.array([R_E]), np.array([0.0])); xr, yr = section_xy(np.array([R_E]), np.array([40.0]))
    ax.plot(xs, ys + 45, marker='*', ms=10, mfc='#e66101', mec='k', mew=0.5, ls='none', clip_on=False, zorder=5)
    ax.plot(xr, yr + 45, marker='v', ms=7, mfc='k', mec='k', ls='none', clip_on=False, zorder=5)
    ax.text(0.01, 0.99, lab, transform=ax.transAxes, va='top', fontsize=9)
cax = fig.add_axes(inch(0.30 * FW, Y_CB1, 0.40 * FW, CBH)); cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cb.set_label('$K_{v_P}$ (10$^{-16}$ s m$^{-3}$)', labelpad=1)
cb.ax.tick_params(labelsize=7, length=2); cb.outline.set_linewidth(0.5)
ax = axs[2]; k, r_km, lon = K['3d']; X, Y = section_xy(r_km, lon); s_ = sector(r_km, lon)
im2 = ax.pcolormesh(X[s_], Y[s_], 100 * dvp3[s_], cmap='RdBu', norm=TwoSlopeNorm(0, -30, 30), shading='gouraud', rasterized=True)
frame(ax); ax.text(0.01, 0.99, '(c)', transform=ax.transAxes, va='top', fontsize=9)
cax2 = fig.add_axes(inch(0.10 * FW, Y_CB2, 0.30 * FW, CBH)); cb2 = fig.colorbar(im2, cax=cax2, orientation='horizontal', extend='both')
cb2.set_label('$\\delta v_P$ (%)', labelpad=1); cb2.ax.tick_params(labelsize=7, length=2); cb2.outline.set_linewidth(0.5)
# (d) test: traveltime change against the size of the v_p perturbation; lines = kernel prediction (linear), symbols = direct runs
ax = fig.add_axes(inch(0.62 * FW, 0.45, 0.36 * FW, Y_C + SECH - 0.45))
amp = np.array([-1.0, -0.5, 0.5, 1.0]); xx = np.linspace(-1.15, 1.15, 50)
for dim, mk, col, lab in (('1d', 'o', 'k', '1-D'), ('3d', 's', '#c1272d', '3-D')):
    z = np.load(f'{D}/test_kernels_{dim}.npz'); cases = [str(c) for c in z['cases']]; meas = z['meas']; P = float(z['pred_blob'])
    ys_ = [1e3 * meas[cases.index(f'{dim}_{c}')][0] for c in ('B_m100', 'B_m050', 'B_p050', 'B_p100')]
    ax.plot(xx, 1e3 * P * xx, color=col, lw=0.8)
    ax.plot(amp, ys_, marker=mk, ms=4.5, mfc='white' if dim == '1d' else col, mec=col, mew=0.8, ls='none', label=lab)
ax.axhline(0, color='0.6', lw=0.5, zorder=0); ax.axvline(0, color='0.6', lw=0.5, zorder=0)
ax.set_xlabel('$\\delta v_P$ at the center (%)'); ax.set_ylabel('$\\delta T$ of P (ms)')
ax.legend(frameon=False, fontsize=7, loc='upper right', handletextpad=0.2); ax.text(-0.28, 1.0, '(d)', transform=ax.transAxes, va='top', fontsize=9)
for s_ in ('top', 'right'): ax.spines[s_].set_visible(False)
ax.tick_params(length=2.5, width=0.5)
fig.savefig(OUT + '.png', dpi=600, bbox_inches='tight', pad_inches=0.02); fig.savefig(OUT + '.pdf', bbox_inches='tight', pad_inches=0.02)
from PIL import Image
_im = Image.open(OUT + '.png'); _a = _im.getchannel('A') if _im.mode == 'RGBA' else None
assert _a is None or _a.getextrema() == (255, 255), 'non-opaque pixels in the PNG'
_im = _im.convert('RGB'); _im.save(OUT + '.png', dpi=(600, 600)); _im.save(OUT + '.jpg', quality=95, subsampling=0, dpi=(600, 600))
print('wrote', OUT + '.png/.pdf/.jpg (RGB, 600 dpi)')
