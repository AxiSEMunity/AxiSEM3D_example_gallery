"""Fig. 4a data from the sphere run (Nr = 360): the phi = 0 slice of U_z (same extraction as paper_figure_tools/package.sh,
ex04) and the y = 0 slices of the input v_p and v_s perturbations.
usage: package_ex4a.py RUN_OUTPUT_ELEMENTS_DIR MODEL_NC OUTDIR"""
import os, sys
import numpy as np, h5py
SRC, NC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
fs = sorted(f for f in os.listdir(SRC) if 'axisem3d_synthetics.nc' in f)
coords, waves, time, nag = [], [], None, None
for f in fs:
    with h5py.File(os.path.join(SRC, f), 'r') as d:
        if nag is None:
            nag = int(d['list_na_grid'][0]); time = d['data_time'][:]
        coords.append(d['list_element_coords'][:]); waves.append(d[f'data_wave__NaG={nag}'][:])
coords = np.concatenate(coords); wave = np.concatenate(waves)
nelem, ngll = coords.shape[0], coords.shape[1]
coords_sz = coords.reshape(nelem * ngll, 2).astype('float32')
wzt = np.empty((nelem * ngll, len(time)), dtype='float32')
for ie in range(nelem):                                       # slice 0 = phi 0, channel 2 = U_z
    wzt[ie * ngll:(ie + 1) * ngll, :] = wave[ie, 0, :, 2, :]
np.savez_compressed(f'{OUT}/cartesian_cross_section.npz', coords_sz=coords_sz, wave_zt=wzt, times=time.astype('float32'))
print(f'cartesian slice: {len(fs)} files, NaG {nag}, {coords_sz.shape[0]} points, {len(time)} samples ({time[0]:.2f}..{time[-1]:.2f} s), max |U_z| {np.abs(wzt).max():.3e}')
with h5py.File(NC, 'r') as d:
    x, y, z = d['x'][()], d['y'][()], d['z'][()]; iy = int(np.argmin(np.abs(y)))
    np.savez_compressed(f'{OUT}/cartesian_model.npz', x=x.astype('float32'), z=z.astype('float32'), yslice=np.array([float(y[iy])]),
                        dvp_xz=d['vp'][:, iy, :].astype('float32'), dvs_xz=d['vs'][:, iy, :].astype('float32'), drho_xz=d['rho'][:, iy, :].astype('float32'))
    print(f'model slice y = {y[iy]:.1f} m: dvp {d["vp"][:, iy, :].min():.2f}, dvs {d["vs"][:, iy, :].min():.2f}, drho {d["rho"][:, iy, :].min():.2f}')
