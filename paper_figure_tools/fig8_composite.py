#!/usr/bin/env python3
"""Fig. 8 of the AxiSEM3D release paper, as one Matplotlib figure: (a) Mars (Example 8), (b) Enceladus (Example 9),
(c) the Sun (Example 10).

Mode v6 makes the figure in the paper: every quantity of rows (a) and (b) is low-passed at the mesh period
(fig8_lowpass.py: Mars 20 s, Enceladus 2 s); row (c) is from the Sun run with a 382-s source half duration, whose
records are low-passed at the 381.6-s mesh period. Modes v3repro, v4 and v5 are earlier versions of the figure and
need data that are not distributed.

Data (folder FIG8_ROOT, default: the current folder)
  mars/       run of 08_atmosphere_Mars_global/axisem3d_mars_atm/input (records and element group mars_xs_fluid)
  enceladus/  run of 09_icesheet_regional/input
  sun/        run of 10_helioseismology_sun/input
  sun_xs.npz  Sun snapshots, written by dump_sun_v6.py from sun/output/elements/sun_cross_section

Usage: python fig8_composite.py v6 out.pdf [--cache DIR]
Afterwards run pdf_plain_images.py on the output (Matplotlib writes palette images, which some PDF conversions
show as black boxes). The paper used a Helvetica subset; without that font file the script uses DejaVu Sans.
"""
import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import font_manager
from matplotlib.ticker import MultipleLocator
from scipy.interpolate import griddata

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)        # fig8_lowpass.py (python -I does not add the script folder)
ROOT = os.environ.get('FIG8_ROOT', '.')   # folder with mars/, enceladus/, sun/ and sun_xs.npz (see the docstring)
EX = BF = ROOT + '/submitted_figure_data'   # modes v3repro and v4 only (data of the submitted figure; not distributed)
MARS_RUN = ROOT + '/mars'
ENC_COMPLETE = ROOT + '/enceladus'
SUN_STATIONS = ROOT + '/sun/output/stations/stationsforplotting'
FONT = os.path.join(HERE, 'fonts', 'Helvetica_fig8_subset.ttf')
FONT_CHARS = set(' ()-./0123456789=AEMRTUabcdefghiklmnoprstuxz|°')
DEG = '°'

SOURCES = {
    'v3repro': dict(
        mars_st=MARS_RUN + '/output/stations/MARS',               # 27 Jun run, stopped at 274 s
        mars_info=BF + '/ex07_mars/data/mars.txt',
        enc_st=EX + '/08_icesheet_regional/output/stations/JW_stations',   # 27 Jun run, stopped at 59 s
        enc_info=BF + '/ex08_icesheet/data/3Dstations.txt'),
    'v4': dict(
        mars_st=MARS_RUN + '/output/stations/MARS',               # 1800 s
        mars_info=MARS_RUN + '/input/mars.txt',
        enc_st=ENC_COMPLETE + '/output/stations/JW_stations',     # 1000 s
        enc_info=ENC_COMPLETE + '/input/3Dstations.txt'),
}

# v5: rerun of job 16807512 (Bouchet job 33571414) with the atmosphere cross-section written every 0.486 s
MARS_XS_DENSE = MARS_RUN + '/output/elements/mars_xs_fluid'
SOURCES['v5'] = dict(SOURCES['v4'], mars_xs=MARS_XS_DENSE)

# v6: as v5; row (c) from the Sun rerun (Bouchet job 33590392): the production run of job 16719245 with the source
# half duration raised from 1 s to 382 s (= mesh period 381.6 s), so no content below the mesh period; same snapshot
# times as v4/v5; snapshots dumped by dump_sun_v6.py
SUN_HD382 = ROOT + '/sun'
SOURCES['v6'] = dict(SOURCES['v5'], sun_npz=ROOT + '/sun_xs.npz',
                     sun_st=SUN_HD382 + '/output/stations/stationsforplotting',
                     sun_lp=(381.6, 800.0), sun_xmax=11200.0)   # records low-passed at the mesh period (end taper 800 s; axis to 11,200 s)

# v5 and v6 (low-passed at the mesh period): choices
V5_MARS_SNAP_T = [150.0, 300.0, 600.0, 900.0, 1200.0, 1500.0]  # s; the acoustic front at 0.5, 1.1, 2.3, 3.4, 4.6, 5.7 deg
V5_MARS_XS_COMMON_SCALE = False    # each cross-section scaled by its own 99th percentile (as in v3/v4)
V5_ENC_MAP_T = [31.0, 46.0, 62.0, 78.0, 93.0, 109.0]          # s; main arrival at the receiver rings 15-52.5 deg
V5_ENC_RS_TMAX = 200.0             # s; end of the time axis of the Enceladus record section

VERMILLION = '#D55E00'      # Okabe-Ito "secondary" of pub_style.py (source stars)

STYLE = {   # pub_style.STYLE of the Bouchet project (the parts that matter here) + TrueType Helvetica
    'font.family': 'sans-serif', 'font.sans-serif': ['Helvetica'],
    'font.size': 10, 'axes.linewidth': 0.8,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4,
    'lines.linewidth': 1.0, 'axes.grid': False, 'axes.spines.top': False, 'axes.spines.right': False,
    'figure.dpi': 300, 'savefig.dpi': 300,            # rasterized layers at 300 dpi, as submitted
    'pdf.fonttype': 42, 'ps.fonttype': 42,             # TrueType text
    'axes.unicode_minus': False, 'mathtext.default': 'regular',
    'figure.constrained_layout.use': False,
    'image.composite_image': True,
}


def setup_font():
    if not os.path.exists(FONT):   # the paper used a Helvetica subset, which is not distributed
        STYLE["font.sans-serif"] = ["DejaVu Sans"]; mpl.rcParams.update(STYLE); return
    font_manager.fontManager.addfont(FONT)
    mpl.rcParams.update(STYLE)
    found = font_manager.findfont(font_manager.FontProperties(family='Helvetica'), fallback_to_default=False)
    assert os.path.samefile(found, FONT), found


def check_text(fig):
    bad = set()
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible():
            bad |= set(t.get_text()) - FONT_CHARS
    if bad:
        sys.exit(f'characters not in the Helvetica subset: {sorted(bad)!r}')


# ---------------------------------------------------------------- data
def read_stations(folder, cache=None):
    """AxiSEM3D ASCII station folder -> time, {name: (nt, 3) array} (columns U1 U2 U3 = R T Z)."""
    fn = None
    if cache:
        os.makedirs(cache, exist_ok=True)
        fn = os.path.join(cache, folder.strip('/').replace('/', '__') + '.npz')
        if os.path.exists(fn):
            d = np.load(fn)
            return d['time'], dict(zip(d['names'], d['data']))
    time = np.loadtxt(os.path.join(folder, 'data_time.ascii'))
    out = {}
    for f in sorted(glob.glob(os.path.join(folder, '*.ascii'))):
        nm = os.path.basename(f)[:-len('.ascii')]
        if nm == 'data_time':
            continue
        a = np.loadtxt(f)
        assert a.shape == (time.size, 3), (f, a.shape)
        out[nm] = a
    if fn:
        np.savez(fn, time=time, names=np.array(list(out)), data=np.array(list(out.values())))
    return time, out


def epicentral_distance_deg(src_lat, src_lon, st_lat, st_lon):
    la1, lo1, la2, lo2 = map(np.radians, (src_lat, src_lon, st_lat, st_lon))
    cosd = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    return float(np.degrees(np.arccos(np.clip(cosd, -1.0, 1.0))))


def record_section(ax, time, recs, gain, lw, scale=None, xmax=None):
    """Vertical-component record section offset by distance. scale None: each trace normalized to its
    own peak (times gain); scale = number: every trace multiplied by the same factor (one amplitude
    scale for the whole section). xmax: right end of the time axis (default: end of the record)."""
    for dist, tr in recs:
        if scale is None:
            y = tr / max(np.max(np.abs(tr)), 1e-30) * gain + dist      # (same arithmetic as v3/v4)
        else:
            y = tr * scale + dist
        ax.plot(time, y, color='#000000', lw=lw)
    dists = np.array([r[0] for r in recs])
    ax.set_xlim(time[0], time[-1] if xmax is None else xmax)
    ax.set_ylim(dists.min() - gain * 1.5, dists.max() + gain * 1.5)
    ax.set_xlabel('Time after source origin (s)', fontsize=16)
    ax.set_ylabel(f'Epicentral distance ({DEG})', fontsize=16)
    ax.tick_params(labelsize=14)


def mars_xs_lowpassed(folder, times, dmax=6.3):
    """v5: atmosphere cross-section of the dense Mars rerun: U3 at every output point (GLL points
    [0, 2, 4] x [0, 2, 4] of each element, as in job 16807512) within dmax deg, low-passed in time at
    the 20-s mesh period, taken at the sample nearest to each of `times`.
    Returns delta (deg), altitude (km), W (points, len(times)), the times shown, and the time axis."""
    import h5py
    from fig8_lowpass import lowpass, MESH_PERIOD, END_TAPER
    R_MARS = 3.3895e6
    fs = sorted(glob.glob(os.path.join(folder, 'axisem3d_synthetics.nc.rank*')))
    if not fs:
        sys.exit(f'no element output in {folder}')
    S, Z, W, t = [], [], [], None
    for f in fs:
        with h5py.File(f, 'r') as h:
            tt = h['data_time'][:]
            if t is None:
                t = tt
            assert np.array_equal(t, tt), f
            co = h['list_element_coords'][:]
            if co.size == 0:
                continue
            wv = [k for k in h.keys() if k.startswith('data_wave')][0]
            w = h[wv][:]
            ne, ng = co.shape[0], co.shape[1]
            S.append(co[:, :, 0].ravel()); Z.append(co[:, :, 1].ravel())
            W.append(w[:, 0, :, 0, :].reshape(ne * ng, -1))
    s = np.concatenate(S); z = np.concatenate(Z); W = np.concatenate(W, axis=0)
    delta = np.degrees(np.arctan2(s, z)); alt = (np.sqrt(s ** 2 + z ** 2) - R_MARS) / 1e3
    sel = delta <= dmax
    dts = np.diff(t)
    assert dts.max() - dts.min() < 1e-3 * dts.mean(), 'element output not evenly sampled'
    Wf = lowpass(W[sel], float(dts.mean()), MESH_PERIOD['mars'], END_TAPER['mars'])
    idx = [int(np.argmin(np.abs(t - x))) for x in times]
    return delta[sel].astype('f4'), alt[sel].astype('f4'), Wf[:, idx], t[idx], t


def time_ticks(ax, time, step):
    """Ticks every `step` s, but none at the right end of the axis (its label would leave the page)."""
    ax.set_xticks(np.arange(0, time[-1] - 1, step))


# ---------------------------------------------------------------- panel (a): Mars
def panel_a(fig, cell, src, mode, cache, log):
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=cell, width_ratios=[1.25, 1.0], wspace=0.22)
    if mode in ('v5', 'v6'):
        # left: the dense rerun, low-passed at 20 s (time series of every point, then the snapshot)
        delta, alt, W, tsnap, txs = mars_xs_lowpassed(src['mars_xs'], V5_MARS_SNAP_T)
        log(f'(a) Mars cross-sections, low-passed at 20 s: {W.shape[0]} points, element output every '
            f'{np.diff(txs).mean():.3f} s, {txs.size} samples, t = {txs[0]:.2f} ... {txs[-1]:.2f} s')
        p99 = np.array([np.percentile(np.abs(W[:, j]), 99.0) for j in range(6)])
        log('(a) 99th percentile of |Uz| per snapshot: ' + ', '.join(f'{x:.2e}' for x in p99))
        vcommon = float(p99.max())
    else:
        # left: six unrolled atmosphere cross-sections (dump of output_xs, job 16807512)
        d = np.load(BF + '/mars_xs_snapshots.npz')
        delta, alt, W, tsnap = d['delta'], d['alt'], d['W'], d['t_snap']
    g = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=gs[0], hspace=0.3, wspace=0.1)
    axs = []
    for j in range(6):
        ax = fig.add_subplot(g[j // 2, j % 2])
        val = W[:, j]
        vmax = max(float(np.percentile(np.abs(val.astype(np.float64)), 99.0)), 1e-30)
        if mode in ('v5', 'v6') and V5_MARS_XS_COMMON_SCALE:
            vmax = vcommon
        ax.scatter(delta, alt, c=val, s=2.5, cmap='RdBu_r', vmin=-vmax, vmax=vmax, rasterized=True,
                   linewidths=0, marker='s')
        ax.plot(0, 0, '*', ms=11, color=VERMILLION, mec='k', mew=0.6, zorder=10, clip_on=False)
        ax.set_xlim(0, 6); ax.set_ylim(0, 60)
        ax.set_xticks([0, 2, 4, 6]); ax.set_yticks([0, 30, 60])
        ax.tick_params(labelsize=14)
        if mode in ('v5', 'v6'):      # the filtered wavefront now reaches the old label position: title above
            ax.set_title(f't = {tsnap[j]:.0f} s', fontsize=15, pad=4)
        else:
            ax.text(0.98, 0.86, f't = {tsnap[j]:.0f} s', transform=ax.transAxes, ha='right', va='top', fontsize=15)
        if j < 4:
            ax.set_xticklabels([])
        if j % 2:
            ax.set_yticklabels([])
        if j == 2:
            ax.set_ylabel('Altitude (km)', fontsize=16)
        axs.append(ax)
    log(f'(a) Mars cross-sections: t = {", ".join(f"{t:.1f}" for t in tsnap)} s '
        + ('(dense rerun, low-passed)' if mode in ('v5', 'v6') else '(dump of output_xs)'))
    # common x label of the six cross-sections, centered under the left block
    pt = 1.0 / (fig.get_figheight() * 72.0)
    fig.text(0.5 * (axs[4].get_position().x0 + axs[5].get_position().x1), axs[4].get_position().y0 - 44.936 * pt,
             f'Epicentral distance ({DEG})', ha='center', va='baseline', fontsize=16)
    # color bar of the snapshots (inside the first panel; a figure axes drawn after the six panels)
    p = axs[0].get_position()
    cax = fig.add_axes([p.x0 + 0.5 * p.width, p.y0 + 0.22 * p.height, 0.46 * p.width, 0.09 * p.height])
    cb = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-1, 1), cmap='RdBu_r'),
                      cax=cax, orientation='horizontal')
    cb.set_ticks([-1, 0, 1])
    cax.tick_params(labelsize=10, length=2, pad=1.5)
    cax.set_title('Uz / max|Uz|', fontsize=12, pad=3)
    # right: surface record section, stations 0-5 deg
    time, data = read_stations(src['mars_st'], cache)
    info = np.loadtxt(src['mars_info'], dtype=str, skiprows=3)
    name2ll = {f'{r[1]}.{r[0]}': (float(r[2]), float(r[3])) for r in info}
    recs = sorted(((epicentral_distance_deg(0., 0., *name2ll[nm]), a[:, 2]) for nm, a in data.items()
                   if nm in name2ll), key=lambda r: r[0])
    recs = [r for r in recs if r[0] <= 5.0]
    if mode in ('v5', 'v6'):                 # low-passed at the 20-s mesh period; each trace to its own peak
        from fig8_lowpass import lowpass, MESH_PERIOD, END_TAPER
        dt = float(np.mean(np.diff(time)))
        recs = [(d, lowpass(tr, dt, MESH_PERIOD['mars'], END_TAPER['mars'])) for d, tr in recs]
    dists = np.array([r[0] for r in recs])
    gain = max(dists.max() - dists.min(), 1.0) / 18.0
    axR = fig.add_subplot(gs[1])
    record_section(axR, time, recs, gain, lw=0.7)
    if mode == 'v3repro':
        axR.xaxis.set_major_locator(MultipleLocator(60))
    else:
        time_ticks(axR, time, 300)
    axR.set_title('Record section', fontsize=16)
    log(f'(a) Mars record section: {len(recs)} stations {dists.min():.0f}-{dists.max():.0f} deg, '
        f'{time.size} samples, t = {time[0]:.2f} ... {time[-1]:.2f} s ({src["mars_st"]})')


# ---------------------------------------------------------------- panel (b): Enceladus
def enceladus_snapshot_times(time, uz):
    """Rule of the submitted figure: six times evenly spaced from the first arrival at the nearest
    receiver to the first arrival at the farthest receiver plus 10 % of that span (arrival = first
    sample above 5 % of the receiver's peak), capped at the end of the record. With the 59-s record
    the cap cut the series at 59 s; with the complete record it ends when the wave has crossed the
    whole network."""
    def onset(tr, frac=0.05):
        m = np.max(np.abs(tr))
        return time[int(np.argmax(np.abs(tr) > frac * m))] if m > 0 else time[-1]
    arr = np.array([onset(tr) for tr in uz])
    t0 = float(arr.min())
    t1 = float(min(arr.max() + 0.10 * (arr.max() - arr.min()), time[-1]))
    return np.linspace(t0, t1, 6), arr


def panel_b(fig, cell, src, mode, cache, log):
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=cell, width_ratios=[1.15, 1.0], wspace=0.22)
    time, data = read_stations(src['enc_st'], cache)
    info = np.loadtxt(src['enc_info'], dtype=str, skiprows=2)
    name2ll = {f'{r[1]}.{r[0]}': (float(r[2]), float(r[3])) for r in info}
    # receivers sorted by distance (record section); ties keep the station-file order
    records = sorted((epicentral_distance_deg(0., 0., *name2ll[nm]), nm, a) for nm, a in data.items()
                     if nm in name2ll)
    lon = np.array([name2ll[r[1]][1] for r in records])
    lat = np.array([name2ll[r[1]][0] for r in records])
    uz = np.array([r[2][:, 2] for r in records])                                 # (receiver, time)
    if mode in ('v5', 'v6'):
        # low-passed at the 2-s mesh period; maps: each map scaled by its largest receiver amplitude
        # (one scale per map, no scaling per receiver); record section: one scale for all traces
        from fig8_lowpass import lowpass, MESH_PERIOD, END_TAPER
        uz = lowpass(uz, float(np.mean(np.diff(time))), MESH_PERIOD['enceladus'], END_TAPER['enceladus'])
        its = [int(np.argmin(np.abs(time - x))) for x in V5_ENC_MAP_T]
        uz_n = np.zeros_like(uz)
        for x in its:
            uz_n[:, x] = uz[:, x] / max(np.abs(uz[:, x]).max(), 1e-30)
        log('(b) largest receiver amplitude per map: '
            + ', '.join(f'{time[x]:.0f} s {np.abs(uz[:, x]).max():.2e} m' for x in its))
    else:
        uz_n = uz / np.maximum(np.max(np.abs(uz), axis=1, keepdims=True), 1e-30)  # each receiver / own peak
    # maps: receivers in station-file order (as submitted; the Delaunay triangulation of the regular
    # receiver grid, and so the interpolated field, depends on the order)
    file_names = [f'{r[1]}.{r[0]}' for r in info if f'{r[1]}.{r[0]}' in data]
    pos = {r[1]: i for i, r in enumerate(records)}
    perm = np.array([pos[n] for n in file_names])
    pts = np.column_stack([lon[perm], lat[perm]])
    glon = np.linspace(lon.min() - 3, lon.max() + 3, 240)
    glat = np.linspace(lat.min() - 3, lat.max() + 3, 240)
    GLON, GLAT = np.meshgrid(glon, glat)
    if mode in ('v5', 'v6'):
        six_t = np.array(V5_ENC_MAP_T)
    else:
        six_t, arr = enceladus_snapshot_times(time, uz)
        log(f'(b) Enceladus: {len(records)} receivers; first arrivals {arr.min():.2f}-{arr.max():.2f} s; '
            f'snapshot targets {", ".join(f"{t:.2f}" for t in six_t)} s')
    g = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[0], wspace=0.12, hspace=0.32)
    shown = []
    for k, t in enumerate(six_t):
        ax = fig.add_subplot(g[k // 3, k % 3])
        it = int(np.argmin(np.abs(time - t)))
        field = np.ma.masked_invalid(griddata(pts, uz_n[perm, it], (GLON, GLAT), method='linear'))
        ax.pcolormesh(glon, glat, field, cmap='RdBu_r', vmin=-0.9, vmax=0.9, rasterized=True, shading='auto')
        ax.scatter(pts[:, 0], pts[:, 1], s=6, c='k', zorder=8, linewidths=0)
        ax.plot(0, 0, '*', ms=13, color=VERMILLION, mec='k', mew=0.6, zorder=10)
        ax.set_aspect('equal')
        ax.set_xticks(np.arange(-30, 61, 30)); ax.set_yticks(np.arange(-60, 61, 30))
        ax.tick_params(labelsize=14)
        ax.grid(True, lw=0.5, color='0.8')      # submitted: gray 0.6 at 50 % opacity = gray 0.8 on white
        if k % 3 == 0:
            ax.set_ylabel(f'lat ({DEG})', fontsize=16)
        else:
            ax.set_yticklabels([])
        if k // 3 == 1:
            ax.set_xlabel(f'lon ({DEG})', fontsize=16)
        else:
            ax.set_xticklabels([])
        ax.set_title(f't = {time[it]:.0f} s', fontsize=15)
        shown.append(time[it])
    log(f'(b) Enceladus maps at t = {", ".join(f"{t:.2f}" for t in shown)} s')
    axR = fig.add_subplot(gs[1])
    recs = [(r[0], r[2][:, 2]) for r in records] if mode not in ('v5', 'v6') else [(r[0], uz[k]) for k, r in enumerate(records)]
    dists = np.array([r[0] for r in recs])
    gain = max(dists.max() - dists.min(), 1.0) / 22.0
    if mode in ('v5', 'v6'):
        # one amplitude scale for all traces: the 90th percentile of the receivers' peaks gets the
        # excursion that every trace had in v4 (gain); weaker receivers stay weaker
        p90 = float(np.percentile(np.abs(uz).max(axis=1), 90.0))
        record_section(axR, time, recs, gain, lw=0.45, scale=gain / p90, xmax=V5_ENC_RS_TMAX)
        time_ticks(axR, np.array([time[0], V5_ENC_RS_TMAX]), 50)
        log(f'(b) record section: one scale, 90th percentile of the receiver peaks = {p90:.2e} m -> {gain:.2f} deg')
    else:
        record_section(axR, time, recs, gain, lw=0.45)
    if mode == 'v3repro':
        axR.xaxis.set_major_locator(MultipleLocator(15))
    elif mode == 'v4':
        time_ticks(axR, time, 200)
    axR.yaxis.set_major_locator(MultipleLocator(15))
    log(f'(b) Enceladus record section: {len(recs)} receivers {dists.min():.1f}-{dists.max():.1f} deg, '
        f'{time.size} samples, t = {time[0]:.2f} ... {time[-1]:.2f} s ({src["enc_st"]})')


# ---------------------------------------------------------------- panel (c): the Sun (unchanged)
def panel_c(fig, cell, log, src=None):
    src = src or {}
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=cell, width_ratios=[1.15, 1.0], wspace=0.22)
    d = np.load(src.get('sun_npz', BF + '/sun_xs_snapshots.npz'))
    S, Z, W, tsnap = d['S'], d['Z'], d['W'], d['t_snap']
    Sfull = np.concatenate([S, -S]) / 1e8
    Zfull = np.concatenate([Z, Z]) / 1e8
    R_SUN = 6.96e8
    theta = np.linspace(0, 2 * np.pi, 400)
    g = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=gs[0], wspace=0.05, hspace=0.22)
    for j in range(6):
        ax = fig.add_subplot(g[j // 3, j % 3])
        val = W[:, j]
        vmax = max(float(np.percentile(np.abs(val.astype(np.float64)), 98.0)), 1e-30)
        ax.scatter(Sfull, Zfull, c=np.concatenate([val, val]), s=0.5, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   rasterized=True, linewidths=0, marker='s')
        ax.plot(R_SUN / 1e8 * np.cos(theta), R_SUN / 1e8 * np.sin(theta), color='#555555', lw=0.7)
        ax.plot(0, (R_SUN - 5e7) / 1e8, '*', ms=12, color=VERMILLION, mec='k', mew=0.6, zorder=10)
        ax.set_xlim(-7.4, 7.4); ax.set_ylim(-7.4, 7.4)
        ax.set_aspect('equal'); ax.axis('off')
        ax.set_title(f't = {tsnap[j] / 3600:.1f} h', fontsize=15)
        if j == 0:
            ax.plot([-1.0, 1.0], [-5.9, -5.9], color='k', lw=2.2, solid_capstyle='butt')
            ax.text(0, -5.2, '200 Mm', ha='center', va='bottom', fontsize=12)
    sun_st = src.get('sun_st', SUN_STATIONS)
    time_st = np.loadtxt(os.path.join(sun_st, 'data_time.ascii'))
    DS = max(1, len(time_st) // 3000)
    ts = time_st[::DS]
    recs = []
    for f in sorted(glob.glob(os.path.join(sun_st, 'a0.d*.ascii'))):
        dist = int(os.path.basename(f).split('.d')[1].split('.ascii')[0])
        x = np.loadtxt(f, usecols=[2])
        if src.get('sun_lp'):
            from fig8_lowpass import lowpass
            x = lowpass(x - x[0], float(time_st[1] - time_st[0]), *src['sun_lp'])
        recs.append((dist, x[::DS]))
    recs.sort(key=lambda r: r[0])
    dists = np.array([r[0] for r in recs])
    gain = (dists.max() - dists.min()) / 30.0
    axR = fig.add_subplot(gs[1])
    record_section(axR, ts, recs, gain, lw=0.5, xmax=src.get('sun_xmax'))
    axR.xaxis.set_major_locator(MultipleLocator(2500))
    axR.yaxis.set_major_locator(MultipleLocator(40))
    log(f'(c) Sun: snapshots t = {", ".join(f"{t / 3600:.2f}" for t in tsnap)} h; record section '
        f'{len(recs)} stations, {ts.size} samples (every {DS}th), t = {ts[0]:.2f} ... {ts[-1]:.2f} s')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('mode', choices=sorted(SOURCES))
    ap.add_argument('out')
    ap.add_argument('--cache', default=None, help='folder for .npz copies of the station files (optional)')
    a = ap.parse_args()
    if os.path.exists(a.out):
        sys.exit(f'{a.out} exists; not overwritten')
    setup_font()
    log = print
    fig = plt.figure(figsize=(11, 15.4))
    fig.set_layout_engine('none')
    rows = gridspec.GridSpec(3, 1, figure=fig, left=0.075, right=0.985, top=0.975, bottom=0.035, hspace=0.36)
    panel_a(fig, rows[0], SOURCES[a.mode], a.mode, a.cache, log)
    panel_b(fig, rows[1], SOURCES[a.mode], a.mode, a.cache, log)
    panel_c(fig, rows[2], log, SOURCES[a.mode])
    pt = 1.0 / (fig.get_figheight() * 72.0)
    for k, letter in enumerate('abc'):           # panel letters, 9.029 pt above the top of each row
        fig.text(0.012, rows[k].get_position(fig).y1 + 9.029 * pt, f'({letter})', fontsize=22,
                 ha='left', va='baseline')
    check_text(fig)
    fig.savefig(a.out)
    log(f'wrote {a.out} (Matplotlib {mpl.__version__})')


if __name__ == '__main__':
    main()
