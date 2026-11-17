"""Adjoint source of the kernel example (Example 7, Fig. 7).

Measurement: cross-correlation traveltime of the vertical displacement at the receiver (0 N, 40 E) in the window
425-475 s with 15-s cosine tapers on both sides (410-490 s):
    dT = - int w v0 ds dt / int w v0^2 dt        (dT > 0 = later; v0 = d/dt s0; s0 = reference trace)
Adjoint force (vertical, at the receiver), as in the published example (tau = 2 s, so kernels give tau * dT):
    F(t_adj) = tau * w(t) * v0(t) / int w v0^2 dt,   t = TSTAR - t_adj.
Timing, exact by construction:
  - every solver time is an integer multiple of dt (main.cpp: step0 = floor(t0 / dt)); TSTAR = N * dt with N = round(600 / dt);
    the force is given at t_adj = TSTAR - t_k for the forward steps t_k, so it is sampled exactly on the adjoint steps;
  - one extra zero sample half a step before the first one fixes the first adjoint step at N - k_last - 1 regardless of
    rounding; k_last is chosen so that the element outputs of both runs (every 10 steps) pair one-to-one
    (k_last = n_f0 - 1 mod 10, n_f0 = absolute step of the first forward output);
  - the force starts and ends with zeros, so FIRST_LAST padding pads with zeros;
  - the adjoint record (663 s) reaches beyond the start of the forward run (TSTAR - t_first_forward = 662.7 s).
Reference trace: the receiver output (every step) of the clean forward run redo/fwd<dim>.
usage: make_adjoint_source.py     (writes redo/bwd<dim>/input/STF/UZ.txt and redo/adjoint_source_<dim>.npz)
"""
import glob, os
import numpy as np, h5py
W = os.environ.get('KERNEL_WORK', os.getcwd())   # folder with redo/, redo_full/ and runs/ (see README)
T1, T2, TAPER, TAU, STRIDE = 425.0, 475.0, 15.0, 2.0, 10
trapz = getattr(np, 'trapezoid', None) or np.trapz


def read_rec(run):
    for fn in sorted(glob.glob(f'{run}/output/stations/Receiver_DT/axisem3d_synthetics.nc.rank*')):
        with h5py.File(fn, 'r') as f:
            st = [b''.join(s).decode('utf-8', 'ignore').strip('\x00 ') for s in f['list_station'][:]]
            if not st: continue
            ch = [b''.join(c).decode('utf-8', 'ignore').strip('\x00 ') for c in f['list_channel'][:]]
            t = f['data_time'][:].astype(float); d = f['data_wave'][:].astype(float)
        return t, d[0, ch.index('U3'), :]
    raise FileNotFoundError(run)


def window(t):
    w = ((t >= T1) & (t <= T2)).astype(float)
    a = (t >= T1 - TAPER) & (t < T1); w[a] = 0.5 * (1 - np.cos(np.pi * (t[a] - (T1 - TAPER)) / TAPER))
    b = (t > T2) & (t <= T2 + TAPER); w[b] = 0.5 * (1 + np.cos(np.pi * (t[b] - T2) / TAPER))
    return w


for dim in ('3d', '1d'):
    t, s0 = read_rec(f'{W}/redo/fwd{dim}'); dt = t[1] - t[0]
    n_abs = np.rint(t / dt).astype(int); assert np.allclose(n_abs * dt, t, atol=1e-6 * dt + 1e-9), 'forward times off the dt grid'
    N = int(round(600.0 / dt)); TSTAR = N * dt
    v0 = np.gradient(s0, dt); w = window(t); Nn = trapz(w * v0 ** 2, dx=dt); g = TAU * w * v0 / Nn
    nz = np.where(w > 0)[0]; k0, k1 = nz[0] - 1, nz[-1] + 1          # one zero sample on each side
    n_f0 = n_abs[0]                                                  # first forward output step (outputs every STRIDE steps from here)
    while (n_abs[k1] - (n_f0 - 1)) % STRIDE: k1 += 1                 # k_last = n_f0 - 1 (mod STRIDE)
    k = np.arange(k0, k1 + 1)
    tadj = (TSTAR - t[k])[::-1]; F = g[k][::-1]
    tadj = np.r_[tadj[0] - 0.5 * dt, tadj]; F = np.r_[0.0, F]
    assert F[0] == 0 and F[1] == 0 and F[-1] == 0 and np.all(np.diff(tadj) > 0)
    step0_adj = int(np.floor(tadj[0] / dt)); assert step0_adj == N - n_abs[k1] - 1
    assert (N - n_f0 - step0_adj) % STRIDE == 0, 'element outputs would not pair'
    out = f'{W}/redo/bwd{dim}/input/STF'; os.makedirs(out, exist_ok=True)
    np.savetxt(f'{out}/UZ.txt', np.c_[tadj, F], fmt='%.10f', delimiter='\t')
    np.savez(f'{W}/redo/adjoint_source_{dim}.npz', t=t, s0=s0, v0=v0, w=w, N=Nn, g=g, tadj=tadj, F=F, TSTAR=TSTAR, TAU=TAU, dt=dt)
    print(f'{dim}: dt {dt:.4f}; TSTAR = {N} dt = {TSTAR:.4f} s; window weights > 0 {t[nz[0]]:.3f}..{t[nz[-1]]:.3f} s; int w v0^2 {Nn:.6e}; '
          f'force: {len(F)} samples, adjoint time {tadj[0]:.4f}..{tadj[-1]:.4f} s, peak {np.abs(F).max():.2f}; '
          f'adjoint run starts at step {step0_adj} = {step0_adj*dt:.4f} s; forward outputs start at step {n_f0}; '
          f'adjoint record must reach {TSTAR - t[0]:.3f} s')
