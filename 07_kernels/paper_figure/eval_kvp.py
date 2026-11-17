"""v_p kernel of Example 7 (Fig. 7) on a set of points.

K_lambda(x) = int div s(x, t) * div s~(x, t) dt   (s = forward field, s~ = adjoint field reversed in time, t = TSTAR - t_adj)
The relative v_p kernel shown in the new Fig. 7 and used for the test:
    K_lnvp(x) = 2 rho vp^2 K_lambda / tau          dT = int K_lnvp dln(vp) dV        (s m^-3; dT > 0 = later)
(v_p changes at fixed v_s and rho change lambda only: d lambda = 2 rho vp d vp.) div s is the trace of the gradient, a
scalar, so the two runs (each in the RTZ frame of its own axis) need no rotation. Stored here: K_lambda and the 1-D rho,
vp at the points; the 3-D vp is added by the integration script (S362ANI x 5).
Forward and adjoint element outputs are paired exactly (make_adjoint_source.py aligns them; checked here).
usage: eval_kvp.py --dim 3d|1d --grid blob|blob2|slice --part I --nparts N --out FILE
"""
import argparse, os, time
import numpy as np

W = os.environ.get('KERNEL_WORK', os.getcwd())   # folder with redo/, redo_full/ and runs/ (see README)


def grid_points(name):
    """(r [m], lat [deg], lon [deg]) flattened, and the grid shape"""
    if name == 'blob':      # the test region: blob center (0 N, 20 E, 960 km depth) +- 3 sigma; 0.5 deg x 25 km
        r = np.arange(4661.0e3, 6161.0e3 + 1, 25.0e3); la = np.arange(-15.0, 15.0001, 0.5); lo = np.arange(5.0, 35.0001, 0.5)
    elif name == 'blob2':   # the same region at half the resolution (convergence check)
        r = np.arange(4661.0e3, 6161.0e3 + 1, 50.0e3); la = np.arange(-15.0, 15.0001, 1.0); lo = np.arange(5.0, 35.0001, 1.0)
    elif name == 'slice':   # source-receiver plane (equator), CMB to surface, -15 to 55 deg
        r = np.arange(3485.0e3, 6371.0e3 + 1, 5.0e3); la = np.array([1e-4]); lo = np.arange(-15.0, 55.0001, 0.1)
    else:
        raise ValueError(name)
    R, LA, LO = np.meshgrid(r, la, lo, indexing='ij')
    return R.ravel(), LA.ravel(), LO.ravel(), R.shape


def integ(a, dt):
    return 0.5 * dt * (a[..., 1:] + a[..., :-1]).sum(axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dim', required=True); ap.add_argument('--grid', required=True)
    ap.add_argument('--part', type=int, default=0); ap.add_argument('--nparts', type=int, default=1)
    ap.add_argument('--out', required=True); ap.add_argument('--batch', type=int, default=800)
    ap.add_argument('--root', default='redo', help='folder under W with fwd<dim>, bwd<dim> (redo or redo_full)')
    ap.add_argument('--full', action='store_true', help='element output with all 25 GLL points (axk_full reader)')
    a = ap.parse_args()
    import sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # axk_full.py
    if a.full:
        from axk_full import ElementOutputFull as EO
    else:
        from axikernels.core.handlers.element_output import ElementOutput as EO
    fwd = EO(path_to_element_output=f'{W}/{a.root}/fwd{a.dim}/output/elements')
    bwd = EO(path_to_element_output=f'{W}/{a.root}/bwd{a.dim}/output/elements')
    adj = np.load(f'{W}/redo/adjoint_source_{a.dim}.npz')
    TSTAR, TAU = float(adj['TSTAR']), float(adj['TAU'])
    tf = np.asarray(fwd.element_groups_info[fwd.element_groups[0]]['metadata']['data_time'], float)
    tb = np.asarray(bwd.element_groups_info[bwd.element_groups[0]]['metadata']['data_time'], float)
    dt = tf[1] - tf[0]; assert abs((tb[1] - tb[0]) - dt) < 1e-6
    tr = (TSTAR - tb)[::-1]                                  # forward time of the reversed adjoint samples
    k0 = int(np.ceil((max(tf[0], tr[0]) - tf[0]) / dt - 1e-6)); k1 = int(np.floor((min(tf[-1], tr[-1]) - tf[0]) / dt + 1e-6))
    tc = tf[k0:k1 + 1]; j0 = int(round((tc[0] - tr[0]) / dt)); off = (tc[0] - tr[0]) / dt - j0
    print(f'forward {tf[0]:.4f}..{tf[-1]:.4f}, adjoint reversed {tr[0]:.4f}..{tr[-1]:.4f}, common {tc[0]:.4f}..{tc[-1]:.4f} '
          f'({tc.size} samples), dt {dt:.4f}, pairing offset {off:+.2e} samples', flush=True)
    assert abs(off) < 1e-3, 'outputs not paired: rebuild the adjoint source'
    assert tc[0] <= tf[0] + 1e-6, 'the adjoint run does not reach the start of the forward run'
    R, LA, LO, shape = grid_points(a.grid)
    idx = np.arange(R.size)[a.part::a.nparts]
    pts = np.column_stack([R[idx], np.radians(LA[idx]), np.radians(LO[idx])])
    n = len(idx); K = np.full(n, np.nan); t00 = time.time()
    for b0 in range(0, n, a.batch):
        sl = slice(b0, min(b0 + a.batch, n))
        df = np.nan_to_num(fwd.load_data(pts[sl], channels=['GRR', 'GTT', 'GZZ'], in_deg=False))[:, :, k0:k1 + 1]
        db = np.nan_to_num(bwd.load_data(pts[sl], channels=['GRR', 'GTT', 'GZZ'], in_deg=False))[:, :, ::-1][:, :, j0:j0 + tc.size]
        K[sl] = integ(df.sum(axis=1) * db.sum(axis=1), dt)
        print(f'  {sl.stop} / {n} points ({time.time() - t00:.0f} s)', flush=True)
    rho = fwd._load_material_property(pts, 'rho'); vp = fwd._load_material_property(pts, 'vp')
    np.savez(a.out, grid=a.grid, dim=a.dim, shape=np.array(shape), index=idx, r=R[idx], lat=LA[idx], lon=LO[idx],
             K_lambda=K, rho1d=rho, vp1d=vp, TSTAR=TSTAR, TAU=TAU, tc0=tc[0], tc1=tc[-1])
    print('saved', a.out, 'K_lambda range', np.nanmin(K), np.nanmax(K), 'NaN', int(np.isnan(K).sum()))


if __name__ == '__main__':
    main()
