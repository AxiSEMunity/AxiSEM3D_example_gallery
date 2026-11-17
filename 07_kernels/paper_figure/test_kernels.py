"""Kernel test of Example 7 (Fig. 7d): kernel prediction against direct simulation.

Measurement (the one the adjoint source encodes): dT = - int w v0 (s - s0) dt / int w v0^2 dt, vertical displacement at the
receiver, window 425-475 s with 15-s cosine tapers; s0 = reference run (ref1d / ref3d, same setup as the kernel runs).
Also: the adjoint functional M = -(1/tau) int F(TSTAR - t) (s - s0) dt with the force F actually used (must equal dT),
and the least-squares time shift (nonlinear check).
(1) v_p test: Gaussian v_p perturbation (1 % at the center; runs at +-0.5 % and +-1 %);
    prediction = int K_lnvp dln(vp) dV,  K_lnvp = 2 rho vp^2 K_lambda / tau  (vp = local model value: PREM, or PREM x
    (1 + 5 dvs_S362ANI) in 3-D; rho from PREM in both), on the 0.5-deg x 25-km grid (and the 1-deg x 50-km grid as a check).
usage: test_kernels.py DIM          (DIM = 3d or 1d)
"""
import glob, os, sys
import numpy as np, h5py
W = os.environ.get('KERNEL_WORK', os.getcwd())   # folder with redo/, redo_full/ and runs/ (see README)
dim = sys.argv[1]; ROOT = sys.argv[2] if len(sys.argv) > 2 else 'redo'      # redo or redo_full (kernel files)
trapz = getattr(np, 'trapezoid', None) or np.trapz
adj = np.load(f'{W}/redo/adjoint_source_{dim}.npz'); TSTAR, TAU = float(adj['TSTAR']), float(adj['TAU'])
CLAT, CLON, SH, RC, SR = 0.0, 20.0, 5.0, 6371.0e3 - 960.0e3, 250.0e3


def angdist(lat1, lon1, lat2, lon2):
    a1, a2 = np.radians(lat1), np.radians(lat2); dl = np.radians(lon2 - lon1)
    h = np.sin((a2 - a1) / 2) ** 2 + np.cos(a1) * np.cos(a2) * np.sin(dl / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1))))


def read_rec(run):
    for fn in sorted(glob.glob(f'{run}/output/stations/Receiver_DT/axisem3d_synthetics.nc.rank*')):
        with h5py.File(fn, 'r') as f:
            st = [b''.join(s).decode('utf-8', 'ignore').strip('\x00 ') for s in f['list_station'][:]]
            if not st: continue
            ch = [b''.join(c).decode('utf-8', 'ignore').strip('\x00 ') for c in f['list_channel'][:]]
            return f['data_time'][:].astype(float), f['data_wave'][0, ch.index('U3'), :].astype(float)
    return None


# ---------------------------------------------------------------- direct runs
t, s0 = read_rec(f'{W}/runs/ref{dim}'); dt = t[1] - t[0]
tf, sf = read_rec(f'{W}/{ROOT}/fwd{dim}')
assert np.allclose(t, tf), 'time axes of ref and fwd differ'
print(f'reference trace: ref{dim} (Nr 50) against the kernel forward run fwd{dim}: max |diff| / max |s0| = {np.abs(sf - s0).max() / np.abs(s0).max():.2e}')
w = np.interp(t, adj['t'], adj['w']); v0 = np.gradient(s0, dt); N0 = trapz(w * v0 ** 2, dx=dt)
F = np.interp(TSTAR - t, adj['tadj'], adj['F'], left=0.0, right=0.0)
from scipy.interpolate import CubicSpline
cs = CubicSpline(t, s0); wm = w > 0
meas = {}
print(f'\n{"case":12s} {"dT (ms)":>11s} {"M adj (ms)":>11s} {"LS shift (ms)":>14s} {"max|ds|/max|s0|":>16s}')
for c in sorted(os.path.basename(os.path.dirname(p)) for p in glob.glob(f'{W}/runs/{dim}_*/output')):
    d = read_rec(f'{W}/runs/{c}')
    if d is None or len(d[1]) != len(s0): print(f'{c:12s} incomplete'); continue
    ds = d[1] - s0
    dT = -trapz(w * v0 * ds, dx=dt) / N0; M = -trapz(F * ds, dx=dt) / TAU
    tw = t[wm]; sh = np.linspace(-1.0, 1.0, 2001); err = [np.sum(w[wm] * (d[1][wm] - cs(tw - x)) ** 2) for x in sh]
    i = int(np.argmin(err)); sh2 = np.linspace(sh[max(i - 1, 0)], sh[min(i + 1, len(sh) - 1)], 2001)
    ls = sh2[int(np.argmin([np.sum(w[wm] * (d[1][wm] - cs(tw - x)) ** 2) for x in sh2]))]
    meas[c] = (dT, M, ls)
    print(f'{c:12s} {1e3*dT:+11.4f} {1e3*M:+11.4f} {1e3*ls:+14.4f} {np.abs(ds).max()/np.abs(s0).max():16.3e}')

# ---------------------------------------------------------------- (1) v_p test
vp3 = None
if dim == '3d':
    import netCDF4
    from scipy.interpolate import RegularGridInterpolator
    with netCDF4.Dataset(f'{W}/redo/fwd3d/input/S362ANI_radius.nc') as ds_:
        rr = np.array(ds_['radius'][:], float); la = np.array(ds_['latitude'][:], float); lo = np.array(ds_['longitude'][:], float)
        dvs = np.array(ds_['dvs'][:], float)
    print(f'\nS362ANI_radius.nc: dvs {dvs.shape}, radius {rr.min():.0f}..{rr.max():.0f}, lat {la.min()}..{la.max()}, lon {lo.min()}..{lo.max()}')
    o = np.argsort(rr); vp3 = RegularGridInterpolator((rr[o], la, lo), dvs[o], bounds_error=False, fill_value=0.0)
pred = {}
for g in ('blob', 'blob2'):
    fns = sorted(glob.glob(f'{W}/{ROOT}/kernels/kvp_{dim}_{g}_part*.npz'))
    if not fns: print('no kernel files for', g); continue
    z0 = np.load(fns[0]); shape = tuple(int(x) for x in z0['shape']); n = int(np.prod(shape))
    A = {k: np.full(n, np.nan) for k in ('K_lambda', 'rho1d', 'vp1d', 'r', 'lat', 'lon')}
    for fn in fns:
        z = np.load(fn)
        for k in A: A[k][z['index']] = z[k]
    miss = int(np.isnan(A['K_lambda']).sum())
    vp = A['vp1d'] if vp3 is None else A['vp1d'] * (1 + 5.0 * vp3(np.c_[A['r'], A['lat'], A['lon']]))
    Kln = 2 * A['rho1d'] * vp ** 2 * A['K_lambda'] / TAU
    dv = 0.01 * np.exp(-0.5 * (angdist(A['lat'], A['lon'], CLAT, CLON) / SH) ** 2) * np.exp(-0.5 * ((A['r'] - RC) / SR) ** 2)
    r_ = A['r'].reshape(shape)[:, 0, 0]; la_ = A['lat'].reshape(shape)[0, :, 0]; lo_ = A['lon'].reshape(shape)[0, 0, :]
    wt = lambda x: np.r_[0.5, np.ones(len(x) - 2), 0.5] * (x[1] - x[0])
    W3 = (wt(r_) * r_ ** 2)[:, None, None] * (wt(np.radians(la_)) * np.cos(np.radians(la_)))[None, :, None] * wt(np.radians(lo_))[None, None, :]
    P = float(np.nansum((Kln * dv).reshape(shape) * W3)); pred[g] = P
    vol_dv = float(np.sum(dv.reshape(shape) * W3)) / 0.01
    print(f'\nv_p kernel on grid {g} {shape}: {len(fns)} parts, missing {miss}; K_lnvp {np.nanmin(Kln):.3e}..{np.nanmax(Kln):.3e} s/m^3; '
          f'vp/vp1d {np.nanmin(vp/A["vp1d"]):.3f}..{np.nanmax(vp/A["vp1d"]):.3f}; blob volume (integral of dv/1%) {vol_dv:.3e} m^3')
    print(f'   predicted dT for +1 %: {1e3*P:+.4f} ms')
print(f'\n{"v_p test":12s} {"measured dT":>12s} {"predicted":>12s} {"pred/meas":>10s} {"(coarse grid)":>14s}')
for c, f_ in (('B_p050', 0.5), ('B_m050', -0.5), ('B_p100', 1.0), ('B_m100', -1.0)):
    k = f'{dim}_{c}'
    if k in meas and 'blob' in pred:
        p = f_ * pred['blob']; pc = f_ * pred.get('blob2', np.nan)
        print(f'{k:12s} {1e3*meas[k][0]:+12.4f} {1e3*p:+12.4f} {p/meas[k][0]:10.4f} {1e3*pc:+14.4f}')
for a_ in ('050', '100'):
    p, m = f'{dim}_B_p{a_}', f'{dim}_B_m{a_}'
    if p in meas and m in meas:
        print(f'   +-{int(a_)/100:.1f} %: odd part {0.5e3*(meas[p][0]-meas[m][0]):+.4f} ms, even part {0.5e3*(meas[p][0]+meas[m][0]):+.4f} ms')
np.savez(f'{W}/{ROOT}/test_kernels_{dim}.npz', cases=np.array(sorted(meas)), meas=np.array([meas[c] for c in sorted(meas)]),
         pred_blob=pred.get('blob', np.nan), pred_blob2=pred.get('blob2', np.nan))
