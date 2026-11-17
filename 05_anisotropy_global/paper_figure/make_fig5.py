#!/usr/bin/env python3
"""Fig. 5 (Example 5, SKS splitting, 2012-07-03 Cook Strait event). (a) Radial and transverse records of the
anisotropic and isotropic 5-s runs, zero-phase low-pass filtered at the 5-s mesh period; (b) map of fast directions and
delay times measured with SplitRacer on the anisotropic synthetics, with the observations of Long et al. (2009).
Run outputs are read from the two run folders of this example (EX05_ROOT, default: the example folder).

usage: python make_fig5.py <synthetic splitting file (lon lat phi dt station, 2 header lines)>
                              <output stem> [--fullpage]
       --fullpage: keep the submitted page size instead of tight edges (layout check only)
"""
import os, sys, glob
HERE = os.path.dirname(os.path.abspath(__file__))
# needs cartopy, shapely, pyproj (pip install cartopy)
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as pe
import matplotlib.ticker as mticker
import cartopy
# Natural Earth shapefiles: cartopy downloads them on first use
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.geodesic import Geodesic
from obspy.taup import TauPyModel

# ------------------------------------------------------------------ inputs
EX = (os.environ.get('EX05_ROOT', os.path.join(HERE, '..')) + '/'   # folder of this example
      '2012-07-03_paper_example_50s/')
OLIV_DIR = EX + 'sim_US32_olivineE_5s/output/stations/JW_stations'     # job 16616637
ISO_DIR = EX + 'sim_US32_isotropic_5s/output/stations/JW_stations'     # job 16688018
STATIONS = EX + 'sim_US32_isotropic_5s/input/eq3.txt'
REAL_FILE = os.path.join(HERE, 'data', 'splitting', 'Long_splitting.txt')  # Long et al. (2009)

# ------------------------------------------------------------------ style (as the submitted figure)
plt.rcParams.update({
    'font.family': 'serif', 'font.serif': ['DejaVu Serif'], 'mathtext.fontset': 'dejavuserif',
    'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 9, 'figure.dpi': 300, 'savefig.dpi': 300, 'axes.linewidth': 0.8,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'lines.linewidth': 1.0,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.grid': False,
    'axes.spines.top': False, 'axes.spines.right': False, 'figure.constrained_layout.use': False,
})
COLORS = {'primary': '#0072B2', 'secondary': '#D55E00', 'quaternary': '#E69F00', 'gray': '#636363'}
C_TRACE = '#5D9CE3'      # anisotropic traces of panel (a) in the submitted figure
C_ISO = 'black'
C_REAL, C_SYN = COLORS['secondary'], COLORS['primary']
LAND, COAST, BORDER = '#efe6cf', COLORS['gray'], '#3a3a3a'
FS = dict(title=16, axlabel=15, tick=13, legend=12.5, panel=18)

# page of the submitted figure and positions measured in it (pt, origin bottom-left)
PW, PH = 906.377, 616.548
RECT_R = (57.309, 156.645, 172.8, 380.16)
RECT_T = (235.869, 156.645, 172.8, 380.16)
RECT_MAP = (521.53, 119.34, 374.62, 475.2)
RECT_CARD = (464.97, 116.06, 123.03, 134.48)          # inset card axes (patch drawn at 0.02-0.98)
RECT_GLOBE = (471.19, 127.35, 109.95, 109.95)

# map
LON_MIN, LON_MAX, LAT_MIN, LAT_MAX = -122.3, -115.3, 40.8, 47.2
LAT0 = 0.5 * (LAT_MIN + LAT_MAX)
KM_PER_S = 48.0
LW = 1.5
SOURCE_LAT, SOURCE_LON = -40.023, 173.756

# panel (a)
EVENT = (-40.023, 173.756)
DEPTH = 229.8
WIN_PRE = 15.0
SHIFT = 5.0
MODEL = TauPyModel('prem')


def frac(r):
    x, y, w, h = r
    return [x / PW, y / PH, w / PW, h / PH]


# ================================================================== panel (a)  (make_composite_figure.py)
def gc_deg(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    return np.degrees(np.arccos(np.clip(
        np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo2 - lo1), -1, 1)))


def load_distances():
    n2d = {}
    for line in open(STATIONS):
        c = line.split()
        if len(c) >= 4 and not line.startswith('#'):
            d = gc_deg(*EVENT, float(c[2]), float(c[3]))
            n2d[f'{c[1]}.{c[0]}'] = d
            n2d[c[0]] = d
    return n2d


def phase_time(d, phase):
    arr = MODEL.get_travel_times(source_depth_in_km=DEPTH, distance_in_degree=d, phase_list=[phase])
    return arr[0].time if arr else None


MESH_PERIOD_S, END_TAPER_S = 5.0, 25.0     # Example 5 mesh period (minimum_period 4.98 s); end taper as in the audit


def lowpass_cols(x, dt):
    """Zero-phase Butterworth low-pass of order 4 at the mesh period (as for Fig. 8, paper_figure_tools/fig8_lowpass.py),
    applied to each column (time along axis 0) after a half-cosine taper over the last END_TAPER_S seconds."""
    from scipy.signal import butter, sosfiltfilt
    x = np.asarray(x, dtype=np.float64); n = int(round(END_TAPER_S / dt)); w = np.ones(x.shape[0])
    w[-n:] = 0.5 * (1.0 + np.cos(np.linspace(0.0, np.pi, n)))
    sos = butter(4, 1.0 / MESH_PERIOD_S, btype='low', fs=1.0 / dt, output='sos')
    return sosfiltfilt(sos, x * w[:, None], axis=0)


def load_scenario(d, n2d):
    time = np.loadtxt(os.path.join(d, 'data_time.ascii'))
    dt = float(time[1] - time[0])
    recs = {}
    for f in sorted(glob.glob(os.path.join(d, '*.ascii'))):
        nm = os.path.basename(f)[:-6]
        if nm == 'data_time' or nm not in n2d:
            continue
        recs[nm] = (n2d[nm], lowpass_cols(np.loadtxt(f), dt))   # panel (a): resolved content only
    return time, recs


def measure_offset(time, recs, half=12.0):
    offs = []
    for nm, (dist, data) in recs.items():
        t_sks = phase_time(dist, 'SKS')
        if t_sks is None:
            continue
        w = (time >= t_sks - half) & (time <= t_sks + half)
        if not w.any():
            continue
        env = np.maximum(np.abs(data[:, 0][w]), np.abs(data[:, 1][w]))
        if env.max() <= 0:
            continue
        offs.append((time[w] - t_sks)[int(np.argmax(env))])
    return float(np.median(offs)) if offs else 0.0


def panel_records(axR, axT):
    n2d = load_distances()
    t_o, rec_o = load_scenario(OLIV_DIR, n2d)
    t_i, rec_i = load_scenario(ISO_DIR, n2d)
    common = sorted(set(rec_o) & set(rec_i), key=lambda nm: rec_o[nm][0])
    if len(common) > 32:
        idx = np.linspace(0, len(common) - 1, 32).astype(int)
        common = [common[i] for i in idx]
    dists = [rec_o[nm][0] for nm in common]
    span = max(max(dists) - min(dists), 1.0)
    gain = span / 22.0
    offset = measure_offset(t_o, rec_o)
    far = max(dists)
    ts_far, tk_far = phase_time(far, 'SKS'), phase_time(far, 'SKKS')
    win_post = 45.0
    if ts_far and tk_far:
        win_post = float(min(max(45.0, (tk_far - ts_far) + 12.0), 165.0))
    skks_pts = []
    for nm in common:
        dist = rec_o[nm][0]
        t_sks = phase_time(dist, 'SKS')
        if t_sks is None:
            continue
        t_align = t_sks + offset
        w = (t_o >= t_align - WIN_PRE - SHIFT) & (t_o <= t_align + win_post)
        if not w.any():
            continue
        tw = t_o[w] - t_align + SHIFT
        Ro = rec_o[nm][1][:, 0][w]
        To = rec_o[nm][1][:, 1][w]
        norm = max(np.max(np.abs(Ro)), 1e-30)
        axR.plot(tw, Ro / norm * gain + dist, color=C_TRACE, lw=0.9, zorder=3)
        wi = (t_i >= t_align - WIN_PRE - SHIFT) & (t_i <= t_align + win_post)
        Ti = rec_i[nm][1][:, 1][wi]
        twi = t_i[wi] - t_align + SHIFT
        axT.plot(twi, Ti / norm * gain + dist, color=C_ISO, lw=1.3, zorder=2)
        axT.plot(tw, To / norm * gain + dist, color=C_TRACE, lw=0.9, zorder=3)
        t_skks = phase_time(dist, 'SKKS')
        if t_skks is not None and -WIN_PRE <= (t_skks - t_sks) <= win_post:
            skks_pts.append((t_skks - t_sks, dist))
    sp = np.array(sorted(skks_pts, key=lambda p: p[1])) if skks_pts else None
    for ax, ttl in ((axR, 'Radial (anisotropic)'), (axT, 'Transverse')):
        ax.axvline(0.0, color=COLORS['gray'], lw=1.0, ls='--', zorder=1)
        if sp is not None:
            ax.plot(sp[:, 0], sp[:, 1], color=COLORS['quaternary'], lw=1.3, ls='--', zorder=1)
        ax.set_xlim(-WIN_PRE, win_post)
        ax.set_ylim(min(dists) - gain * 1.5, max(dists) + gain * 1.5)
        ax.set_title(ttl, fontsize=FS['title'])
        ax.tick_params(labelsize=FS['tick'])
    axT.tick_params(labelleft=False)


# ================================================================== panel (b)  (make_splitting_figure.py)
def load_split(fname, header):
    d = np.genfromtxt(fname, skip_header=header, usecols=(0, 1, 2, 3))
    d = d[~np.isnan(d).any(axis=1)]
    return d[:, 0], d[:, 1], d[:, 2], d[:, 3]


def stick_segments(lon, lat, phi, dt, geod):
    pts = np.column_stack([lon, lat]).astype(float)
    half_m = 0.5 * dt * KM_PER_S * 1000.0
    a = np.asarray(geod.direct(pts, phi, half_m))
    b = np.asarray(geod.direct(pts, phi + 180.0, half_m))
    return [[(b[i, 0], b[i, 1]), (a[i, 0], a[i, 1])] for i in range(len(lon))]


def halo():
    return [pe.withStroke(linewidth=2.2, foreground='white')]


def draw_map(fig, ax, syn_file):
    pc = ccrs.PlateCarree()
    geod = Geodesic()
    rlon, rlat, rphi, rdt = load_split(REAL_FILE, 0)
    slon, slat, sphi, sdt = load_split(syn_file, 2)
    ax.set_extent([LON_MIN, LON_MAX, LAT_MIN, LAT_MAX], crs=pc)
    ax.add_feature(cfeature.LAND.with_scale('50m'), facecolor=LAND, zorder=0)
    # no lakes (no rivers or lake outlines); the submitted map showed
    # Malheur and Harney Lakes as small white patches
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'), edgecolor=COAST, linewidth=0.7, zorder=6)
    ax.add_feature(cfeature.STATES.with_scale('50m'), edgecolor=BORDER, linewidth=1.3, zorder=6)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), edgecolor=BORDER, linewidth=1.6, zorder=6)
    gl = ax.gridlines(crs=pc, draw_labels=True, linewidth=0.3, color='0.8', linestyle=':', zorder=1)
    gl.top_labels = gl.right_labels = False
    gl.xlocator = mticker.FixedLocator([-122, -120, -118, -116])
    gl.ylocator = mticker.FixedLocator([42, 44, 46])
    # the submitted figure has the 122 W gridline but no label for it
    gl.xformatter = mticker.FuncFormatter(lambda v, pos: '' if abs(v + 122) < 1e-6 else f'{abs(v):.0f}°W')
    gl.xlabel_style = {'size': 13}
    gl.ylabel_style = {'size': 13}
    ax.add_collection(LineCollection(stick_segments(rlon, rlat, rphi, rdt, geod), transform=pc,
                                     colors=C_REAL, linewidths=LW, capstyle='round', zorder=4))
    ax.add_collection(LineCollection(stick_segments(slon, slat, sphi, sdt, geod), transform=pc,
                                     colors=C_SYN, linewidths=LW, capstyle='round', zorder=5))
    # CALIFORNIA (-121.6, 41.25) is left out: in the submitted figure it lies under the inset card
    for name, lo, la in [('WASHINGTON', -118.5, 46.95), ('OREGON', -120.7, 42.9), ('IDAHO', -115.9, 45.4),
                         ('NEVADA', -116.6, 41.3)]:
        ax.text(lo, la, name, transform=pc, fontsize=13, color='#33405e', style='italic',
                ha='center', va='center', zorder=7, path_effects=halo())
    # source -> array globe inset (card and globe at the submitted positions)
    bgax = fig.add_axes(frac(RECT_CARD), zorder=20)
    bgax.add_patch(FancyBboxPatch((0.02, 0.02), 0.96, 0.96, transform=bgax.transAxes,
                                  boxstyle='round,pad=0.0,rounding_size=0.04',
                                  facecolor='white', edgecolor='0.6', linewidth=0.6))
    bgax.set_axis_off()

    def _unit(la, lo):
        la, lo = np.radians(la), np.radians(lo)
        return np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    clat0, clon0 = float(np.mean(rlat)), float(np.mean(rlon))
    mid = _unit(SOURCE_LAT, SOURCE_LON) + _unit(clat0, clon0)
    clon = np.degrees(np.arctan2(mid[1], mid[0]))
    clat = np.degrees(np.arctan2(mid[2], np.hypot(mid[0], mid[1])))
    gax = fig.add_axes(frac(RECT_GLOBE), projection=ccrs.Orthographic(clon, clat), zorder=21)
    gax.set_facecolor('white')
    gax.set_global()
    gax.add_feature(cfeature.LAND.with_scale('110m'), facecolor=LAND, zorder=0)
    gax.add_feature(cfeature.COASTLINE.with_scale('110m'), edgecolor=COAST, linewidth=0.3, zorder=1)
    gax.gridlines(linewidth=0.25, color='0.7', linestyle=':')
    gax.plot([SOURCE_LON, clon0], [SOURCE_LAT, clat0], transform=ccrs.Geodetic(), color='0.35', lw=0.8, zorder=2)
    gax.scatter([clon0], [clat0], transform=pc, marker='s', s=22, c=COLORS['primary'],
                edgecolors='k', linewidths=0.3, zorder=4)
    gax.scatter([SOURCE_LON], [SOURCE_LAT], transform=pc, marker='*', s=150, c=COLORS['secondary'],
                edgecolors='k', linewidths=0.4, zorder=5)
    gax.set_title('source → array', fontsize=11, pad=3)
    gax.spines['geo'].set_linewidth(0.5)


def one_second_pt():
    """On-map length (pt) of a 1-s stick, E-W at the mid-latitude (Mercator x scale)."""
    frac_w = (KM_PER_S / (111.32 * np.cos(np.radians(LAT0)))) / (LON_MAX - LON_MIN)
    return frac_w * RECT_MAP[2]


# ================================================================== figure
def text_pt(fig, x, y, s, size, **kw):
    fig.text(x / PW, y / PH, s, fontsize=size, ha=kw.pop('ha', 'left'), va=kw.pop('va', 'baseline'), **kw)


def line_pt(fig, x0, x1, y, color, lw, cap):
    fig.add_artist(Line2D([x0 / PW, x1 / PW], [y / PH, y / PH], transform=fig.transFigure,
                          color=color, lw=lw, solid_capstyle=cap))


def main(syn_file, stem, fullpage=False):
    fig = plt.figure(figsize=(PW / 72.0, PH / 72.0))
    ax_R = fig.add_axes(frac(RECT_R))
    ax_T = fig.add_axes(frac(RECT_T), sharey=ax_R)
    panel_records(ax_R, ax_T)
    proj = ccrs.Mercator(central_longitude=0.5 * (LON_MIN + LON_MAX), min_latitude=LAT_MIN, max_latitude=LAT_MAX)
    ax_b = fig.add_axes(frac(RECT_MAP), projection=proj)
    draw_map(fig, ax_b, syn_file)
    # labels and legends at the positions of the submitted figure (baseline-left anchors, pt)
    text_pt(fig, 16.9023, 258.7256, 'Epicentral distance (°)', FS['axlabel'], rotation=90, rotation_mode='anchor')
    text_pt(fig, 103.5356, 114.7324, 'Time relative to predicted SKS (s)', FS['axlabel'])
    text_pt(fig, 670.22, 88.20, 'Longitude', FS['axlabel'])
    text_pt(fig, 468.43, 324.99, 'Latitude', FS['axlabel'], rotation=90, rotation_mode='anchor')
    text_pt(fig, 6.7021, 600.8721, '(a)', FS['panel'])
    text_pt(fig, 6.7021 + 24.416 * 18, 600.8721, '(b)', FS['panel'])
    # legend (a): 10 pt higher than in the submitted figure (legends closer to the panels)
    line_pt(fig, 139.4416, 161.9416, 94.0992, C_TRACE, 2.6, 'projecting')
    line_pt(fig, 139.4416, 161.9416, 75.35, C_ISO, 2.6, 'projecting')
    text_pt(fig, 171.94, 89.72, 'Anisotropic (olivine)', FS['legend'])
    text_pt(fig, 171.94, 70.97, 'Isotropic (transverse only)', FS['legend'])
    # legend (b): sample lines = on-map length of a 1-s stick, ending where the submitted lines ended
    L1 = one_second_pt()
    text_pt(fig, 618.68, 64.67, 'Stick length = 1 s splitting', FS['legend'])
    line_pt(fig, 626.41 - L1, 626.41, 44.65, C_REAL, LW + 0.6, 'round')
    line_pt(fig, 623.94 - L1, 623.94, 25.16, C_SYN, LW + 0.6, 'round')
    text_pt(fig, 633.91, 40.28, 'Real data — Long et al. (2009)', FS['legend'])
    text_pt(fig, 631.44, 20.78, 'Synthetic — AxiSEM3D (US32)', FS['legend'])
    print(f'1-s stick on the map: {L1:.2f} pt (legend sample line)')
    # Matplotlib fills every axes background with a "no outline" state of stroke opacity 0 (/CA 0),
    # which a journal converter treats as transparency.  The backgrounds are white on a white page
    # (and under the map: land everywhere, the white card under the globe), so they are not drawn.
    for a in fig.axes:
        a.patch.set_visible(False)
    kw = {} if fullpage else dict(bbox_inches='tight', pad_inches=0.5 / 25.4)
    fig.savefig(stem + '.pdf', **kw)
    fig.savefig(stem + '_rgba.png', dpi=600, facecolor='white', **kw)
    plt.close(fig)
    # opaque 8-bit RGB raster twins (no alpha channel)
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    im = Image.open(stem + '_rgba.png')
    a = np.asarray(im)
    assert a.shape[2] == 4 and a[:, :, 3].min() == 255, 'PNG has transparent pixels'
    rgb = Image.fromarray(np.ascontiguousarray(a[:, :, :3]))
    rgb.save(stem + '.png', dpi=(600, 600))
    rgb.save(stem + '.jpg', 'JPEG', quality=95, subsampling=0, optimize=True, dpi=(600, 600))
    os.remove(stem + '_rgba.png')
    print(f'written {stem}.pdf/.png/.jpg  ({rgb.size[0]} x {rgb.size[1]} px)')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], '--fullpage' in sys.argv)
