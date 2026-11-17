#!/usr/bin/env python3
"""
Mars atmosphere -- the acoustic wave TRAPPED in the 60 km atmosphere, shown by
"unrolling" the thin shell: epicentral distance (deg) along x, altitude above the
solid surface (km) up y. The source is an atmospheric pressure pulse, so the energy
rides the surface/atmosphere rather than diving into the interior -- this shell is
where the wavefield actually lives.

Reads the cross-section re-run atmosphere group:
  output_xs/elements/mars_xs_fluid/   (U3 in the 60 km atmosphere, with (s,z) coords)
plus the surface station data in ./data/ for the record section.
"""
import os, sys, glob
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pub_style import apply_style, save_fig, COLORS, epicentral_distance_deg
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

apply_style()
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['figure.constrained_layout.use'] = False  # let fig.subplots_adjust margins apply
mpl.rcParams['pdf.use14corefonts'] = True    # reference core Helvetica (Illustrator-editable)
mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['axes.unicode_minus'] = False
R_MARS = 3.3895e6
RUN = os.path.join(HERE, '..', 'axisem3d_mars_atm', 'output_xs')


def read_group(group):
    fs = sorted(glob.glob(os.path.join(RUN, 'elements', group,
                                        'axisem3d_synthetics.nc.rank*')))
    d0 = Dataset(fs[0]); d0.set_auto_mask(False)
    t = np.asarray(d0.variables['data_time'][:]); d0.close()
    S, Z, W = [], [], []
    for f in fs:
        d = Dataset(f); d.set_auto_mask(False)
        co = np.asarray(d.variables['list_element_coords'][:])
        wv = [k for k in d.variables if k.startswith('data_wave')][0]
        w = np.asarray(d.variables[wv][:])
        ne, ng = co.shape[0], co.shape[1]
        S.append(co[:, :, 0].ravel()); Z.append(co[:, :, 1].ravel())
        W.append(w[:, 0, :, 0, :].reshape(ne * ng, -1))
        d.close()
    return np.concatenate(S), np.concatenate(Z), np.concatenate(W, axis=0), t


s, z, W, t = read_group('mars_xs_fluid')
r = np.sqrt(s**2 + z**2)
delta = np.degrees(np.arctan2(s, z))          # colatitude = epicentral distance (deg)
alt = (r - R_MARS) / 1e3                       # altitude above solid surface (km)
print(f'{s.size:,} atmosphere points, alt {alt.min():.0f}..{alt.max():.0f} km, '
      f'{len(t)} frames')

# how far has the (slow) acoustic wave reached by the end -> set the distance window
dmax = 6.0                       # cut the view at 6 deg epicentral distance
NS = 6
snap_idx = np.linspace(len(t) // 12, len(t) - 1, NS).astype(int)

# ---- surface record section (same station data as the still figures) ----
ST_DIR = os.path.join(HERE, 'data', 'stations')
st_info = np.loadtxt(os.path.join(HERE, 'data', 'mars.txt'), dtype=str, skiprows=3)
name2ll = {f'{r_[1]}.{r_[0]}': (float(r_[2]), float(r_[3])) for r_ in st_info}
rtime = np.loadtxt(os.path.join(ST_DIR, 'data_time.ascii'))
recs = []
for f in sorted(glob.glob(os.path.join(ST_DIR, '*.ascii'))):
    nm = os.path.basename(f).replace('.ascii', '')
    if 'data_time' in f or nm not in name2ll:
        continue
    recs.append((epicentral_distance_deg(0., 0., *name2ll[nm]), np.loadtxt(f)))
recs = [rr for rr in recs if rr[0] <= dmax]       # match the 6 deg cut
recs.sort(key=lambda rr: rr[0])
dists = np.array([rr[0] for rr in recs])
gain = max(dists.max() - dists.min(), 1.0) / 18.0

# ---- figure: 6 unrolled atmosphere snapshots (left) + record section (right) ----
fig = plt.figure(figsize=(11, 5.8)); fig.set_layout_engine('none')
outer = gridspec.GridSpec(1, 2, width_ratios=[1.25, 1.0], wspace=0.16, figure=fig)
gs = gridspec.GridSpecFromSubplotSpec((NS + 1) // 2, 2, subplot_spec=outer[0], hspace=0.38, wspace=0.10)
sel = delta <= dmax * 1.05
for j in range(NS):
    ax = fig.add_subplot(gs[j // 2, j % 2])
    val = W[sel, snap_idx[j]]
    vmax = max(np.percentile(np.abs(val), 99.0), 1e-30)
    ax.scatter(delta[sel], alt[sel], c=val, s=3, cmap='RdBu_r',
               vmin=-vmax, vmax=vmax, rasterized=True, linewidths=0, marker='s')
    ax.plot(0, 0, '*', ms=10, color=COLORS['secondary'], mec='k', mew=0.6, zorder=10)
    ax.set_xlim(0, dmax); ax.set_ylim(0, 60)
    ax.tick_params(labelsize=12)
    ax.text(0.985, 0.82, f't = {t[snap_idx[j]]:.0f} s', transform=ax.transAxes,
            ha='right', va='top', fontsize=13)
    if j < NS - 1:
        ax.set_xticklabels([])
    if j == NS // 2:
        ax.set_ylabel('Altitude (km)', fontsize=17)
ax.set_xlabel('Epicentral distance (°)', fontsize=17)

axR = fig.add_subplot(outer[1])
for dd, dat in recs:
    tr = dat[:, 2] if dat.ndim > 1 else dat
    nv = max(np.max(np.abs(tr)), 1e-30)
    axR.plot(rtime, tr / nv * gain + dd, color=COLORS['mono'], lw=0.5)
axR.set_xlabel('Time after source origin (s)', fontsize=19)
axR.set_ylabel('Epicentral distance (°)', fontsize=19)
axR.set_xlim(rtime[0], rtime[-1]); axR.tick_params(labelsize=17)
axR.set_ylim(dists.min() - gain * 1.5, dmax + gain)
axR.set_title('Surface record section', fontsize=19)

cax = fig.add_axes([0.07, 0.065, 0.33, 0.014])
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                  cax=cax, orientation='horizontal')
cb.set_label('Uz / max|Uz|  (each snapshot normalised)', fontsize=13)
cb.set_ticks([-1, 0, 1]); cb.ax.tick_params(labelsize=8)
fig.suptitle('Mars atmosphere — trapped acoustic wave (unrolled shell) and record section',
             fontsize=20)
fig.subplots_adjust(left=0.07, right=0.97, top=0.92, bottom=0.20)
save_fig(fig, os.path.join(HERE, 'mars_atm_section'), tight=False)
print('Mars atmosphere section done.')
