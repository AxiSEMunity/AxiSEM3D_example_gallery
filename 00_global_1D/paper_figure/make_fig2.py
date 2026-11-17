#!/usr/bin/env python3
"""
Fig. 2 of the AxiSEM3D release paper: maps of the vertical-displacement wavefield over the
conterminous United States after the 2011 Mw 5.8 Virginia earthquake, for

  (a) Example 1   1D PREM                      run ex1_prem_5s
  (b) Example 2a  S362ANI                      run ex2a_s362ani_5s
  (c) Example 2b  S362ANI + CRUST1.0           run ex2b_crust_s362ani_5s

at t = 200, 450, 700 and 950 s. The defaults (times, 15-s low-pass, color limit +-15 um) give the figure
in the paper.

Runs: the 5-s input folders of Examples 1, 2a and 2b in this repository (Gaussian source time function with half
duration 2.5 s), run with AxiSEM3D v2.1.0.

Input: the station group US_grid_0p1deg of each run (regular 0.1-degree grid, longitude -130 to -63, latitude 21.5 to
52.5, vertical displacement U3, 1-s sampling, one NetCDF-4 file per MPI rank):
  <data>/<run>/output/stations/US_grid_0p1deg/axisem3d_synthetics.nc.rank*

Processing: zero-phase low-pass in time (Butterworth, order 4 forward and backward, corner period --lowpass), so that
the shortest wavelengths left in the maps are resolved by the 0.1-degree grid; cubic interpolation in time to the exact
snapshot times; bicubic resampling to a Lambert conformal conic map (standard parallels 33 and 45 N, central meridian
96 W), in which the wavefronts of the 1D model are close to circles.

Run, for example (about 48 GB of memory)
       python make_fig2.py

Output: <out>.png (600 dpi, 8-bit RGB without alpha), <out>.pdf (TrueType text, no transparency, opaque RGB images)
and <out>.jpg (quality 95, no chroma subsampling; the manuscript includes the figure as .jpg), plus a cache of the
filtered snapshots (cache/; use --fresh after new data arrive).
"""
import argparse
import glob
import io
import os

import h5py
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from scipy.ndimage import map_coordinates
from scipy.signal import butter, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- fixed facts of the runs
GROUP = 'US_grid_0p1deg'
LON0, LAT0, DGRID, NLON, NLAT = -130.0, 21.5, 0.1, 671, 311     # grid of setup/stage_runs.sh
SRC_LAT, SRC_LON = 37.91, -77.93                                 # inparam.source.yaml (VIRGINIA_201108231751A)
RUNS = [                                                         # run folder, row label
    ('ex1_prem_5s', '(a) 1D PREM'),
    ('ex2a_s362ani_5s', '(b) 3D mantle: S362ANI'),
    ('ex2b_crust_s362ani_5s', '(c) 3D mantle and crust: S362ANI + CRUST1.0'),
]

# ---------------------------------------------------------------- map: Lambert conformal conic on a sphere
R_EARTH = 6371.0
FRAME = (-2440.0, 2360.0, -1640.0, 1460.0)       # x0, x1, y0, y1 in km: the conterminous US with a small margin


class Lambert:
    def __init__(self, lat1=33.0, lat2=45.0, lon0=-96.0, lat0=39.0):
        p1, p2 = np.radians([lat1, lat2])
        self.n = np.log(np.cos(p1) / np.cos(p2)) / np.log(np.tan(np.pi / 4 + p2 / 2) / np.tan(np.pi / 4 + p1 / 2))
        self.f = np.cos(p1) * np.tan(np.pi / 4 + p1 / 2) ** self.n / self.n
        self.lon0 = lon0
        self.rho0 = self._rho(lat0)

    def _rho(self, lat):
        return R_EARTH * self.f / np.tan(np.pi / 4 + np.radians(lat) / 2) ** self.n

    def forward(self, lon, lat):
        rho = self._rho(np.asarray(lat, float))
        th = self.n * np.radians(np.asarray(lon, float) - self.lon0)
        return rho * np.sin(th), self.rho0 - rho * np.cos(th)

    def inverse(self, x, y):
        rho = np.hypot(x, self.rho0 - y)
        th = np.arctan2(x, self.rho0 - y)
        lat = np.degrees(2 * np.arctan((R_EARTH * self.f / rho) ** (1 / self.n)) - np.pi / 2)
        return self.lon0 + np.degrees(th / self.n), lat


PROJ = Lambert()


# ---------------------------------------------------------------- data
def read_grid(stadir, tmax):
    """Vertical displacement on the grid: time (nt,), field (nt, NLAT, NLON) in m, up to time tmax."""
    files = sorted(glob.glob(os.path.join(stadir, 'axisem3d_synthetics.nc.rank*')))
    if not files:
        raise FileNotFoundError('no rank files in ' + stadir)
    time = None
    for fn in files:                                         # common time axis (a run in progress: the written part)
        with h5py.File(fn, 'r') as f:
            t = f['data_time'][:]
        with np.errstate(invalid='ignore', over='ignore'):
            bad = np.where(~(np.diff(t) > 0) | (np.abs(t[1:]) > 1e7))[0]
        t = t[:int(bad[0]) + 1] if len(bad) else t
        if time is not None and not np.allclose(t[:len(time)], time[:len(t)], rtol=0, atol=1e-6):
            raise ValueError('time axes differ between rank files: ' + fn)
        time = t if time is None or len(t) < len(time) else time
    n = min(len(time), int(np.searchsorted(time, tmax)) + 1)
    time = time[:n]
    field = np.full((n, NLAT, NLON), np.nan, np.float32)
    for fn in files:
        with h5py.File(fn, 'r') as f:
            keys = np.ascontiguousarray(f['list_station'][:])
            keys = keys.view('S%d' % keys.shape[1]).ravel()
            w = f['data_wave'][:, 0, :][:, :n]          # whole rows in one read, then the part needed
        jj = np.array([int(k[4:7]) for k in keys])
        ii = np.array([int(k[8:11]) for k in keys])
        field[:, jj, ii] = w.T
    if np.isnan(field[0]).any():
        raise ValueError('%d grid points missing in %s' % (np.isnan(field[0]).sum(), stadir))
    return time, field


def lowpass_snapshots(time, field, tc, times, pad=150.0):
    """Zero-phase low-pass (corner period tc, s) along time, then the field at the exact times (cubic in time)."""
    dt = float(np.median(np.diff(time)))
    if tc > 0:
        npad = int(round(pad / dt))                          # the field is zero before the origin time
        x = np.concatenate([np.zeros((npad,) + field.shape[1:], np.float32), field], axis=0)
        sos = butter(4, (1.0 / tc) / (0.5 / dt), output='sos')
        y = sosfiltfilt(sos, x, axis=0)[npad:].astype(np.float32)
    else:
        y = field
    snaps = []
    for ts in times:
        k = int(np.searchsorted(time, ts))
        idx = np.clip(np.arange(k - 2, k + 2), 0, len(time) - 1)
        tt = time[idx]
        if len(set(idx)) < 4:                                # at an end of the record: nearest sample
            snaps.append(y[int(np.argmin(np.abs(time - ts)))].astype(np.float64))
            continue
        wgt = [np.prod([(ts - tt[m]) / (tt[j] - tt[m]) for m in range(4) if m != j]) for j in range(4)]
        snaps.append(np.tensordot(np.array(wgt), y[idx].astype(np.float64), axes=(0, 0)))
    return np.array(snaps)


def front_amplitude(grid2d, t, pct):
    """Percentile of |u| in the annulus that holds the surface waves at time t (2.6 t to 4.2 t km from the epicenter).
    The whole-map maximum is no guide: the permanent displacement at the epicenter is about 1 cm."""
    lon, lat = np.meshgrid(LON0 + DGRID * np.arange(NLON), LAT0 + DGRID * np.arange(NLAT))
    la1, lo1, la2, lo2 = map(np.radians, (SRC_LAT, SRC_LON, lat, lon))
    cosd = np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1)
    dist = R_EARTH * np.arccos(np.clip(cosd, -1, 1))
    ring = (dist > 2.6 * t) & (dist < 4.2 * t)
    if ring.sum() < 100:                                     # first snapshot: the waves have not left the source
        return 0.0
    return float(np.percentile(np.abs(grid2d[ring]), pct))


def to_map(grid2d, nx, ny):
    """Bicubic resampling of a (NLAT, NLON) field to the map frame, ny x nx pixels, row 0 = south."""
    x0, x1, y0, y1 = FRAME
    xs = x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx
    ys = y0 + (np.arange(ny) + 0.5) * (y1 - y0) / ny
    lon, lat = PROJ.inverse(*np.meshgrid(xs, ys))
    fj, fi = (lat - LAT0) / DGRID, (lon - LON0) / DGRID
    if fj.min() < 0 or fj.max() > NLAT - 1 or fi.min() < 0 or fi.max() > NLON - 1:
        raise ValueError('map frame leaves the grid')
    return map_coordinates(grid2d, [fj, fi], order=3, mode='nearest')


def read_lines(name):
    """Polylines of a GMT multi-segment file, projected; only segments that touch the frame."""
    x0, x1, y0, y1 = FRAME
    segs, cur = [], []

    def close():
        if len(cur) > 1:
            a = np.array(cur)
            x, y = PROJ.forward(a[:, 0], a[:, 1])
            if x.max() > x0 and x.min() < x1 and y.max() > y0 and y.min() < y1:
                segs.append(np.column_stack([x, y]))

    with open(os.path.join(HERE, 'map_lines', name)) as f:
        for line in f:
            if line[0] in '>#':
                close()
                cur = []
            else:
                a = line.split()
                cur.append((float(a[0]), float(a[1])))
    close()
    return segs


# ---------------------------------------------------------------- figure
def diverging_cmap():
    """ColorBrewer RdBu (reversed: blue = down, red = up) with a pure white center and without its two darkest end
    colors (#053061, #67001f), so that clipped crests stay blue and red instead of near-black."""
    c = ['#2166ac', '#4393c3', '#92c5de', '#d1e5f0', '#ffffff', '#fddbc7', '#f4a582', '#d6604d', '#b2182b']
    return LinearSegmentedColormap.from_list('RdBu_white', c, N=511)


def star_outline(radius):
    """Closed outline of a five-pointed star (km); drawn with an opaque white fill and a 0.5-pt black edge."""
    a = np.pi / 2 + np.arange(11) * np.pi / 5
    r = np.where(np.arange(11) % 2 == 0, radius, 0.381966 * radius)
    return np.column_stack([r * np.cos(a), r * np.sin(a)])


def nice_step(v):
    """Largest of 1, 2, 5 x 10^k not above v."""
    e = 10.0 ** np.floor(np.log10(v))
    return max(s * e for s in (1, 2, 5) if s * e <= v * (1 + 1e-9))


def make_figure(snaps, times, vmax, unit, scale, out, dpi=600):
    plt.rcParams.update({                                   # the style of the paper's other figures (pub_style.py)
        'font.family': 'serif', 'font.size': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
        'axes.linewidth': 0.5, 'mathtext.fontset': 'cm', 'savefig.dpi': dpi, 'figure.dpi': dpi,
    })
    cmap = diverging_cmap()
    nrow, ncol = len(snaps), len(times)
    x0, x1, y0, y1 = FRAME
    # layout in inches; map pixels match the output resolution exactly, so no image is resampled
    width, gap, top, row_lab, row_gap, cb_block = 7.0, 0.045, 0.20, 0.17, 0.05, 0.47
    nx = int(np.floor((width - (ncol - 1) * gap) / ncol * dpi))
    ny = int(round(nx * (y1 - y0) / (x1 - x0)))
    pw, ph = nx / dpi, ny / dpi
    height = top + nrow * (row_lab + ph) + (nrow - 1) * row_gap + cb_block
    fig = plt.figure(figsize=(width, height))

    shore = read_lines('shore_l.txt')
    nation = read_lines('border_national_l.txt')
    state = read_lines('border_state_l.txt')
    sx, sy = PROJ.forward(SRC_LON, SRC_LAT)
    star = star_outline(105.0)

    for r in range(nrow):
        ytop = height - top - r * (row_lab + ph + row_gap)
        fig.text(0.0, (ytop - 0.035) / height, RUNS[r][1], ha='left', va='top', fontsize=8)
        for c in range(ncol):
            ax = fig.add_axes([c * (pw + gap) / width, (ytop - row_lab - ph) / height, pw / width, ph / height])
            img = to_map(snaps[r][c] * scale, nx, ny)
            rgb = (cmap(np.clip(0.5 + 0.5 * img / vmax, 0, 1))[:, :, :3] * 255).round().astype(np.uint8)
            ax.patch.set_visible(False)                        # the opaque image is the background
            ax.imshow(rgb, extent=FRAME, origin='lower', interpolation='none', aspect='auto', zorder=1)
            ax.add_collection(LineCollection(state, colors='#9a9a9a', linewidths=0.25, capstyle='round', zorder=2))
            ax.add_collection(LineCollection(nation, colors='#4d4d4d', linewidths=0.35, capstyle='round', zorder=3))
            ax.add_collection(LineCollection(shore, colors='#4d4d4d', linewidths=0.35, capstyle='round', zorder=3))
            ax.fill(sx + star[:, 0], sy + star[:, 1], facecolor='white', edgecolor='k', lw=0.5, joinstyle='miter', zorder=5)
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
            ax.set_xticks([])
            ax.set_yticks([])
            if r == 0:
                fig.text((c * (pw + gap) + pw / 2) / width, (height - 0.03) / height, 't = %g s' % times[c],
                         ha='center', va='top', fontsize=9)
            if r == 0 and c == 0:                              # scale bar in the Gulf of Mexico
                bx, by, bl = 60.0, y0 + 120.0, 1000.0
                ax.plot([bx, bx + bl], [by, by], color='k', lw=0.8, solid_capstyle='butt', zorder=6)
                for xx in (bx, bx + bl):
                    ax.plot([xx, xx], [by - 30, by + 30], color='k', lw=0.6, zorder=6)
                ax.text(bx + bl / 2, by + 60, '1000 km', ha='center', va='bottom', fontsize=7, zorder=6)

    # color bar: an opaque image plus two triangles for the clipped ends
    cw, ch = 2.6, 0.085
    cax = fig.add_axes([(width - cw) / 2 / width, 0.27 / height, cw / width, ch / height])
    cax.patch.set_visible(False)
    grad = (cmap(np.linspace(0, 1, 512))[None, :, :3] * 255).round().astype(np.uint8)
    cax.imshow(grad, extent=(-vmax, vmax, 0, 1), aspect='auto', interpolation='nearest')
    tri = 0.035 * 2 * vmax
    for sgn, col in ((-1, cmap(0.0)), (1, cmap(1.0))):
        cax.add_patch(Polygon([[sgn * vmax, 0], [sgn * (vmax + tri), 0.5], [sgn * vmax, 1]], closed=True,
                              fc=col[:3], ec='k', lw=0.5, clip_on=False))
    step = nice_step(vmax / 2)
    ticks = np.arange(-np.floor(vmax / step), np.floor(vmax / step) + 1) * step
    cax.set_xticks(ticks)
    cax.set_xticklabels([('%g' % t).replace('-', '−') for t in ticks])
    cax.set_yticks([])
    cax.tick_params(axis='x', length=2.5, width=0.5, labelsize=7.5, pad=2)
    cax.set_xlabel('Vertical displacement (%s)' % unit, fontsize=8.5, labelpad=2)

    buf = io.BytesIO()                                     # Agg writes RGBA: store the PNG as 8-bit RGB without alpha
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', pad_inches=0.02)
    buf.seek(0)
    rgba = np.asarray(Image.open(buf))
    if rgba.shape[2] == 4 and (rgba[:, :, 3] != 255).any():
        raise ValueError('the figure has pixels that are not opaque')
    Image.fromarray(np.ascontiguousarray(rgba[:, :, :3])).save(out + '.png', dpi=(dpi, dpi))
    fig.savefig(out + '.pdf', dpi=dpi, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(out + '.jpg', dpi=dpi, bbox_inches='tight', pad_inches=0.02, pil_kwargs={'quality': 95, 'subsampling': 0})
    plt.close(fig)
    print('saved', out + '.png, .pdf and .jpg; map pixels', nx, 'x', ny, '= %.2f km' % ((x1 - x0) / nx))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--data', default=os.path.join(HERE, '..', 'grid_output'))
    ap.add_argument('--out', default=os.path.join(HERE, '..', 'figures', 'Figure_basic_wavefield_snapshots'))
    ap.add_argument('--lowpass', type=float, default=15.0, help='corner period of the low-pass (s); 0 = none')
    ap.add_argument('--times', type=float, nargs='+', default=[200.0, 450.0, 700.0, 950.0])
    ap.add_argument('--vmax', type=float, default=15.0, help='color limit in display units; 0 = automatic')
    ap.add_argument('--unit', default='um', choices=['um', 'mm', 'm'])
    ap.add_argument('--runs', nargs='+', default=None, help='run folders to use (default: the three of the figure)')
    ap.add_argument('--cache', default=os.path.join(HERE, '..', 'cache'))
    ap.add_argument('--dpi', type=int, default=600)
    ap.add_argument('--fresh', action='store_true', help='read the grid output again instead of the cached snapshots')
    args = ap.parse_args()
    if args.runs:
        byname = dict(RUNS)
        RUNS[:] = [(r, byname.get(r, r)) for r in args.runs]
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)

    snaps = []
    for run, _ in RUNS:
        tag = '%s_lp%g_t%s' % (run, args.lowpass, '-'.join('%g' % t for t in args.times))
        cfile = os.path.join(args.cache, tag + '.npz')
        if os.path.exists(cfile) and not args.fresh:
            s = np.load(cfile)['snaps']
        else:
            time, field = read_grid(os.path.join(args.data, run, 'output', 'stations', GROUP), max(args.times) + 150.0)
            print('%s: %d samples to %.1f s, dt = %.5f s, max |u| = %.3e m'
                  % (run, len(time), time[-1], np.median(np.diff(time)), np.abs(field).max()))
            s = lowpass_snapshots(time, field, args.lowpass, args.times)
            np.savez(cfile, snaps=s.astype(np.float32), times=np.array(args.times), lowpass=args.lowpass)
            del field
        snaps.append(s)
        for ts, a in zip(args.times, s):
            print('  %s  t = %5.0f s  surface-wave annulus: 99.5th percentile of |u| = %7.2f um, max = %7.2f um; whole map max = %.1f um'
                  % (run, ts, 1e6 * front_amplitude(a, ts, 99.5), 1e6 * front_amplitude(a, ts, 100.0), 1e6 * np.abs(a).max()))

    unit, scale = {'um': ('µm', 1e6), 'mm': ('mm', 1e3), 'm': ('m', 1.0)}[args.unit]
    vmax = args.vmax
    if vmax <= 0:      # automatic: the surface-wave front of the last snapshot just about reaches full color
        vmax = 0.8 * max(front_amplitude(s[-1], args.times[-1], 99.5) for s in snaps) * scale
        vmax = float('%.2g' % vmax)                     # two significant digits ('%.1g' rounded 13.6 down to 10)
    print('color limit: +-%g %s' % (vmax, args.unit))
    make_figure(snaps, args.times, vmax, unit, scale, args.out, dpi=args.dpi)


if __name__ == '__main__':
    main()
