"""Fig. 4b data from the plume run with a resolved source (half duration 14 s = mesh period): the phi = 0 slice of U_z,
t <= 520 s (the same extraction as paper_figure_tools/package.sh, ex04). usage: package_ex4b.py RUN_OUTPUT_ELEMENTS_DIR OUTDIR"""
import os, sys
import numpy as np, h5py
SRC, OUT = sys.argv[1], sys.argv[2]
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
for ie in range(nelem):
    wzt[ie * ngll:(ie + 1) * ngll, :] = wave[ie, 0, :, 2, :]
keep = time <= 520.0
np.savez_compressed(f'{OUT}/plume_cross_section.npz', coords_sz=coords_sz, wave_zt=wzt[:, keep], times=time[keep].astype('float32'))
print(f'plume slice: {len(fs)} files, {coords_sz.shape[0]} points, {int(keep.sum())} samples to {time[keep][-1]:.1f} s, max |U_z| {np.abs(wzt[:, keep]).max():.3e}')
