#!/bin/bash
# Package a completed AxiSEM3D simulation into a self-contained plots/ subfolder.
#
# Each subfolder ends up with:  pub_style.py, make_figures.py, data/, and the PDF/PNG figures.
# Highest-resolution run for each example is used (5s where available, else native).
#
# Usage:  ./package.sh <example>     where <example> in:
#           ex00 ex01 ex02 ex03 ex04 ex05 ex06 ex07 ex08   (or 'all' for the ready ones)
#
# Safe to re-run: it overwrites data/ and regenerates figures.

set -uo pipefail

PLOTS_DIR="$(cd "$(dirname "$0")" && pwd)"
EX="$(dirname "$PLOTS_DIR")"        # .../examples

load_py() {
    module purge 2>/dev/null || true
    module load StdEnv GCC/13.3.0 OpenMPI/5.0.3-GCC-13.3.0 \
        SciPy-bundle/2024.05-gfbf-2024a matplotlib/3.9.2-gfbf-2024a 2>/dev/null
}

# copytree SRC DST  — copy a directory's ascii station output (recursive)
copytree() {
    local src="$1" dst="$2"
    if [ -d "$src" ]; then
        mkdir -p "$dst"
        cp -r "$src"/. "$dst"/
        return 0
    fi
    echo "    MISSING: $src" >&2
    return 1
}

# copytree_rtz SRC DST — ENZ -> RTZ conversion for the ex05 meridian SKS runs.
# The solver writes 27-column ENZ station files (col0=East, col1=North, col2=Z,
# then 24 gradient/strain/stress columns). The source sits at the north pole, so
# for these meridian receivers radial = local NORTH and transverse = local EAST.
# make_figures.py expects a 3-column RTZ file (R,T,Z), so per station we keep the
# first three ENZ columns and SWAP cols 0<->1: R=North($2), T=East($1), Z=Z($3).
# data_time.ascii is copied verbatim (the time vector matches the new sampling);
# list_channel.info is rewritten to the 3-col RTZ channel list. Idempotent.
copytree_rtz() {
    local src="$1" dst="$2"
    if [ ! -d "$src" ]; then echo "    MISSING: $src" >&2; return 1; fi
    mkdir -p "$dst"
    local f base
    for f in "$src"/*; do
        base=$(basename "$f")
        if [ "$base" = "data_time.ascii" ]; then
            cp "$f" "$dst/$base"
        elif [ "$base" = "list_channel.info" ]; then
            printf 'R\nT\nZ\n' > "$dst/$base"
        elif [[ "$base" == *.ascii ]]; then
            # R=North(col2), T=East(col1), Z=vertical(col3); drop cols 4..27
            awk '{print $2, $1, $3}' "$f" > "$dst/$base"
        fi
    done
    return 0
}

# copytree_disp SRC DST — like copytree but for ENZ [U,G,E,S,R] 27-column station
# files: keep only the first 3 columns (displacement U_E,U_N,U_Z) to cut download
# size ~9x. data_time.ascii and any *.info are copied verbatim. Figures only need
# the displacement channels, so this is lossless for the plots.
copytree_disp() {
    local src="$1" dst="$2"
    if [ ! -d "$src" ]; then echo "    MISSING: $src" >&2; return 1; fi
    mkdir -p "$dst"
    local f base
    for f in "$src"/*; do
        base=$(basename "$f")
        if [ "$base" = "data_time.ascii" ] || [[ "$base" == *.info ]]; then
            cp "$f" "$dst/$base"
        elif [[ "$base" == *.ascii ]]; then
            awk '{print $1, $2, $3}' "$f" > "$dst/$base"
        else
            cp -r "$f" "$dst/$base"
        fi
    done
    return 0
}

run_fig() {
    local dest="$1"
    load_py
    ( cd "$dest" && python3 make_figures.py )
}

# Conda interpreters for bundles that need deps outside the module stack.
#   axisem3d   : cartopy, xarray, netCDF4, matplotlib, pyyaml  (ex10, sun)
#   axikernels : pandas, h5py, tables, netCDF4, matplotlib     (ex11, ex12)
CONDA_AXISEM3D="$(command -v python3)"   # environment with cartopy, xarray, netCDF4
CONDA_AXIKERN="$(command -v python3)"   # environment with pandas, h5py, tables, netCDF4

# run_fig_py PY DEST — run a bundle's make_figures.py with an explicit interpreter.
run_fig_py() {
    local py="$1" dest="$2"
    if [ ! -x "$py" ]; then echo "    MISSING interpreter: $py" >&2; return 1; fi
    ( cd "$dest" && "$py" make_figures.py )
}

# copy_station_subset SRC DST STATION... — copy only data_time.ascii plus the named
# station .ascii files. The GSN group records ~130 global stations but the figures
# plot only one (IU.ANMO), so copying the whole group would waste ~330 MB.
copy_station_subset() {
    local src="$1" dst="$2"; shift 2
    if [ ! -d "$src" ]; then echo "    MISSING: $src" >&2; return 1; fi
    mkdir -p "$dst"
    cp "$src/data_time.ascii" "$dst/" 2>/dev/null || true
    local st
    for st in "$@"; do
        cp "$src/$st.ascii" "$dst/" 2>/dev/null || echo "    (no station $st in $src)" >&2
    done
    return 0
}

case "${1:-}" in
ex00)
    echo "=== ex00 global 1D (5s) ==="
    D="$PLOTS_DIR/ex00_global_1D"; SRC="$EX/00_global_1D"
    rm -rf "$D/data"; mkdir -p "$D/data"
    cp "$SRC/input_5s/GSN.txt" "$SRC/input_5s/US_TA.txt" "$SRC/input_5s/inparam.source.yaml" "$D/data/"
    # only IU.ANMO is plotted (seismogram); the station map uses GSN.txt for locations
    copy_station_subset "$SRC/output_5s/stations/global_seismic_network_GSN" "$D/data/GSN_stations" IU.ANMO
    copytree "$SRC/output_5s/stations/USArray_transportable"       "$D/data/USArray_stations"
    run_fig "$D"
    ;;
ex01)
    echo "=== ex01 S362ANI (5s) ==="
    D="$PLOTS_DIR/ex01_S362ANI"; SRC="$EX/01_S362ANI_EMC_global"
    rm -rf "$D/data"; mkdir -p "$D/data"
    cp "$SRC/input_5s/US_TA.txt" "$D/data/"
    # only IU.ANMO is plotted (1D-vs-3D seismogram)
    copy_station_subset "$SRC/output_5s/stations/global_seismic_network_GSN" "$D/data/3D_GSN_stations" IU.ANMO
    copytree "$SRC/output_5s/stations/USArray_transportable"       "$D/data/3D_USArray_stations"
    # 1D PREM comparison comes from ex00 5s
    copy_station_subset "$EX/00_global_1D/output_5s/stations/global_seismic_network_GSN" "$D/data/1D_GSN_stations" IU.ANMO || \
        echo "    (1D comparison skipped — ex00 5s not ready)"
    run_fig "$D"
    ;;
ex02)
    echo "=== ex02 regional crust (5s) ==="
    D="$PLOTS_DIR/ex02_regional_crust"; SRC="$EX/02_3d_crust_S362ANI_regional"
    rm -rf "$D/data"; mkdir -p "$D/data"
    cp "$SRC/simu_with_3d_crust_5s/input/US_TA.txt" "$D/data/"
    copytree "$SRC/simu_with_1d_crust_5s/output/stations/USArray_transportable" "$D/data/1D_USArray_stations"
    copytree "$SRC/simu_with_3d_crust_5s/output/stations/USArray_transportable" "$D/data/3D_USArray_stations"
    run_fig "$D"
    ;;
ex03)
    echo "=== ex03 salt body (0.2s) ==="
    D="$PLOTS_DIR/ex03_salt_body"; SRC="$EX/03_salt_body_SEG_local"
    rm -rf "$D/data"; mkdir -p "$D/data"
    copytree "$SRC/simu3D/output/elements/Fourier_coefficients_ocean_floor" "$D/data/ocean_floor"
    copytree "$SRC/simu3D/output/elements/orthogonal_azimuthal_slices"      "$D/data/azimuthal_slices"
    # extract compact VP/VS slices of the SEG/EAGE C3 salt body INPUT model
    # (the 1 GB SEG_C3_SOLID.nc is not shipped; we keep two ~0.4 MB slices)
    load_py
    python3 - "$SRC/simu3D/input/SEG_C3_data/SEG_C3_SOLID.nc" "$D/data/salt_model.npz" <<'PYEOF'
import sys, numpy as np, netCDF4
SRC, OUT = sys.argv[1], sys.argv[2]
nc = netCDF4.Dataset(SRC)
x = nc.variables['x'][:].astype('float32'); y = nc.variables['y'][:].astype('float32')
depth = nc.variables['depth'][:].astype('float32')
VP = nc.variables['VP'][:]; VS = nc.variables['VS'][:]   # (depth, x, y), km/s
salt = VP > 4.4                                          # high-VP salt vs slow sediment
yb = int(np.argmax(salt.sum(axis=(0, 1))))              # y-slice through the most salt
db = int(np.argmax(salt.sum(axis=(1, 2))))              # depth of the widest salt
np.savez_compressed(OUT, x=x, y=y, depth=depth,
    ybest=np.array([yb]), dbest=np.array([db]),
    VP_xz=VP[:, :, yb].astype('float32'), VS_xz=VS[:, :, yb].astype('float32'),
    VP_xy=VP[db, :, :].astype('float32'), VS_xy=VS[db, :, :].astype('float32'))
print(f'  salt model: ybest={yb} (y={y[yb]:.0f} m), dbest={db} (depth={depth[db]:.0f} m)')
PYEOF
    run_fig "$D"
    ;;
ex04)
    echo "=== ex04 release paper (0.08s) ==="
    D="$PLOTS_DIR/ex04_simple_3d_shapes"; SRC="$EX/04_simple_3d_shapes"
    rm -rf "$D/data"; mkdir -p "$D/data"
    # plume wavefield cross-section: extract the phi=0, U_z slice (channel 2) from
    # the 2.1 GB element output and keep only the t<=520 s window the figure uses,
    # giving a ~16 MB npz instead of shipping the raw NetCDF.
    PLUME="$SRC/example_single_plume/output/elements/orthogonal_azimuthal_slices"
    if [ -d "$PLUME" ]; then
        load_py
        python3 - "$PLUME" "$D/data/plume_cross_section.npz" <<'PYEOF'
import sys, os, numpy as np, xarray as xr
SRC, OUT = sys.argv[1], sys.argv[2]
fs = sorted([f for f in os.listdir(SRC) if 'axisem3d_synthetics.nc' in f])
d0 = xr.open_dataset(os.path.join(SRC, fs[0]))
nag = int(d0['list_na_grid'].values[0]); time = d0['data_time'].values
coords = xr.concat([xr.open_dataset(os.path.join(SRC, f))['list_element_coords']
                    for f in fs], dim='dim_element').values
wave = xr.concat([xr.open_dataset(os.path.join(SRC, f))[f'data_wave__NaG={nag}']
                  for f in fs], dim=f'dim_element__NaG={nag}').values
nelem, ngll = coords.shape[0], coords.shape[1]
coords_sz = coords.reshape(nelem * ngll, 2).astype('float32')
wzt = np.empty((nelem * ngll, len(time)), dtype='float32')
for ie in range(nelem):                       # slice 0 = phi=0, channel 2 = U_z
    wzt[ie * ngll:ie * ngll + ngll, :] = wave[ie, 0, :, 2, :]
keep = time <= 520.0                           # the figure's snapshots stop at 500 s
np.savez_compressed(OUT, coords_sz=coords_sz, wave_zt=wzt[:, keep],
                    times=time[keep].astype('float32'))
print(f'  plume slice: {coords_sz.shape[0]} pts, {int(keep.sum())} timesteps (<=520 s)')
PYEOF
    else
        echo "    (plume run not found — skipping plume cross-section)"
    fi
    # plume INPUT velocity model: meridional dV slice + head-depth map (compact)
    PLNC="$SRC/example_single_plume/plume_1.nc"
    if [ -f "$PLNC" ]; then
        load_py
        python3 - "$PLNC" "$D/data/plume_model.npz" <<'PYEOF'
import sys, numpy as np, netCDF4
nc = netCDF4.Dataset(sys.argv[1])
lat = nc.variables['lat'][:].astype('float32'); lon = nc.variables['lon'][:].astype('float32')
depth = nc.variables['depth'][:].astype('float32')
vs = nc.variables['vs'][:]                          # (lat,lon,depth) fractional dV
ilon0 = int(np.argmin(np.abs(lon - 0.0)))           # meridional slice lon=0
idep = int(np.argmin(np.abs(depth - 1200e3)))       # plume-head depth map
np.savez_compressed(sys.argv[2], lat=lat, lon=lon, depth=depth,
    dv_xz=vs[:, ilon0, :].astype('float32'), dv_xy=vs[:, :, idep].astype('float32'),
    idep=np.array([idep]), ilon0=np.array([ilon0]))
print(f'  plume model: head depth {depth[idep]/1e3:.0f} km, dV {vs[:,ilon0,:].min():.2f}..0')
PYEOF
    fi
    # cartesian blob wavefield cross-section: same extraction as the plume above,
    # but pointed at the cartesian-blob element output (phi=0 slice, channel 2 =
    # U_z). The record is short (~15 s) so we keep every timestep; make_figures.py
    # picks its own snapshot window. This mirrors plume_cross_section.npz so the
    # two time-lapse figures form a consistent pair.
    CARTEL="$SRC/example_input_cartesian/output/elements/orthogonal_azimuthal_slices"
    if [ -d "$CARTEL" ]; then
        load_py
        python3 - "$CARTEL" "$D/data/cartesian_cross_section.npz" <<'PYEOF'
import sys, os, numpy as np, xarray as xr
SRC, OUT = sys.argv[1], sys.argv[2]
fs = sorted([f for f in os.listdir(SRC) if 'axisem3d_synthetics.nc' in f])
d0 = xr.open_dataset(os.path.join(SRC, fs[0]))
nag = int(d0['list_na_grid'].values[0]); time = d0['data_time'].values
coords = xr.concat([xr.open_dataset(os.path.join(SRC, f))['list_element_coords']
                    for f in fs], dim='dim_element').values
wave = xr.concat([xr.open_dataset(os.path.join(SRC, f))[f'data_wave__NaG={nag}']
                  for f in fs], dim=f'dim_element__NaG={nag}').values
nelem, ngll = coords.shape[0], coords.shape[1]
coords_sz = coords.reshape(nelem * ngll, 2).astype('float32')
wzt = np.empty((nelem * ngll, len(time)), dtype='float32')
for ie in range(nelem):                       # slice 0 = phi=0, channel 2 = U_z
    wzt[ie * ngll:ie * ngll + ngll, :] = wave[ie, 0, :, 2, :]
np.savez_compressed(OUT, coords_sz=coords_sz, wave_zt=wzt,  # short record: keep all t
                    times=time.astype('float32'))
print(f'  cartesian slice: {coords_sz.shape[0]} pts, {len(time)} timesteps '
      f'({time[0]:.2f}..{time[-1]:.2f} s)')
PYEOF
    else
        echo "    (cartesian element run not found — skipping cross-section)"
    fi
    # cartesian blob seismograms (small) + station list (for the distance record
    # section) + a compact y=0 cross-section of the blob velocity anomaly.
    CART=$(find "$SRC/example_input_cartesian/output/stations" -mindepth 1 -maxdepth 1 -type d | head -1)
    [ -n "$CART" ] && copytree "$CART" "$D/data/cartesian_stations"
    cp "$SRC/example_input_cartesian/stn.txt" "$D/data/cartesian_stn.txt" 2>/dev/null || true
    BLOBNC="$SRC/example_input_cartesian/example_sphere_cartesian.nc"
    if [ -f "$BLOBNC" ]; then
        load_py
        python3 - "$BLOBNC" "$D/data/cartesian_blob.npz" <<'PYEOF'
import sys, numpy as np, netCDF4
nc = netCDF4.Dataset(sys.argv[1])
x = nc.variables['x'][:].astype('float32'); y = nc.variables['y'][:]
z = nc.variables['z'][:].astype('float32')
vs = nc.variables['vs'][:]                 # (x, y, z) fractional dVs, REF1D
iy = int(np.argmin(np.abs(y - 0.0)))       # the y=0 plane through source/stations
np.savez_compressed(sys.argv[2], x=x, z=z, dvs_xz=vs[:, iy, :].astype('float32'),
                    yslice=np.array([float(y[iy])]))
print(f'  cartesian blob: y={y[iy]:.0f} m slice, dVs {vs[:,iy,:].min():.2f}..0')
PYEOF
    fi
    # release paper record section
    cp "$SRC/example_release_paper/input/paperstns.txt" "$D/data/" 2>/dev/null || true
    RP=$(find "$SRC/example_release_paper/output/stations" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -1)
    [ -n "$RP" ] && copytree "$RP" "$D/data/release_paper_stations"
    run_fig "$D"
    ;;
ex05)
    echo "=== ex05 anisotropy (5s) — SKS distance profile, RTZ [U] = 3 cols ==="
    D="$PLOTS_DIR/ex05_anisotropy"; SRC="$EX/05_anisotropy_global"
    mkdir -p "$D/data"
    # The three meridian SKS scenarios (sim1 VTI mesh, sim2 iso-mesh+ani overlay,
    # dm lowermost-mantle tilted Cij) were rerun with a corrected vertical-CLVD
    # source on the bouchet mirror; their fresh output lives there (the in-repo
    # examples tree holds only empty run skeletons). The rerun run dirs do NOT use
    # the _5s suffix. Output is 27-col ENZ, converted to 3-col RTZ via copytree_rtz
    # (R=North, T=East, Z=vertical). Bundle dir names keep their established _5s
    # suffix for figure compatibility. US32 (a real earthquake DC) is NOT rerun and
    # is left exactly as shipped.
    RERUN="$EX/../AxiSEM3D_selected_examples_bouchet/05_anisotropy_global"
    PREM="$RERUN/PREM_anisotropy_w_and_wo_full_Cij_50s"
    # distance->station mapping for the SKS plots (lon-0 SK.D## meridian line)
    cp "$PREM/sim1_ani_prem_mesh/input/stations.txt" "$D/data/sks_stations.txt" 2>/dev/null || true
    # ENZ -> RTZ col-swap+truncate for the three meridian scenarios (overwrites the
    # old buggy data in place; idempotent).
    copytree_rtz "$PREM/sim1_ani_prem_mesh/output/stations/JW_stations" \
                 "$D/data/sim1_ani_prem_mesh_5s"             || echo "    (sim1 ani not ready)"
    copytree_rtz "$PREM/sim2_iso_prem_mesh_plus_ani/output/stations/JW_stations" \
                 "$D/data/sim2_iso_prem_mesh_plus_ani_5s"    || echo "    (sim2 iso not ready)"
    copytree_rtz "$RERUN/deep_mantle_anisotropy_full_Cij_50s/sim_lowermost_mantle_ani/output/stations/JW_stations" \
                 "$D/data/sim_lowermost_mantle_ani_5s"       || echo "    (deep mantle not ready)"
    # US32: real 2012-07-03 NZ event at the US TA array (real DC source, NOT rerun).
    # Leave the shipped 3-col RTZ data in place; only refresh if a run is present.
    copytree "$SRC/2012-07-03_paper_example_50s/sim_US32_olivineE_5s/output/stations/JW_stations" \
             "$D/data/sim_US32_olivineE_5s" 2>/dev/null || echo "    (US32 olivine: keeping shipped data)"
    copytree "$SRC/2012-07-03_paper_example_50s/sim_US32_isotropic_5s/output/stations/JW_stations" \
             "$D/data/sim_US32_isotropic_5s" 2>/dev/null || echo "    (US32 isotropic: keeping shipped data)"
    # US32 station file (eq3.txt) so the figure can compute great-circle distances.
    cp "$SRC/2012-07-03_paper_example_50s/sim_US32_olivineE_5s/input/eq3.txt" \
       "$D/data/us32_stations.txt" 2>/dev/null || true
    run_fig "$D"
    ;;
ex06)
    echo "=== ex06 finite source SFBA (2s) ==="
    D="$PLOTS_DIR/ex06_finite_source"; SRC="$EX/06_finite_source_local"
    rm -rf "$D/data"; mkdir -p "$D/data/stations"
    cp "$SRC/input/STATIONS_OUTPUT.txt" "$D/data/"
    SFBA=$(find "$SRC/output/stations" -mindepth 1 -maxdepth 1 -type d 2>/dev/null | head -1)
    if [ -n "$SFBA" ]; then
        cp "$SFBA"/axisem3d_synthetics.nc.* "$D/data/stations/" 2>/dev/null
        cp "$SFBA"/*.ascii "$D/data/stations/" 2>/dev/null || true
    else
        echo "    MISSING: ex06 station output"
    fi
    run_fig "$D"
    ;;
ex07)
    echo "=== ex07 Mars (20s) ==="
    D="$PLOTS_DIR/ex07_mars"; SRC="$EX/07_atmosphere_Mars_global"
    rm -rf "$D/data"; mkdir -p "$D/data"
    cp "$SRC/axisem3d_mars_atm/input/mars.txt" "$D/data/"
    copytree "$SRC/axisem3d_mars_atm/output/stations/MARS" "$D/data/stations"
    run_fig "$D"
    ;;
ex08)
    echo "=== ex08 icesheet (2s) ==="
    D="$PLOTS_DIR/ex08_icesheet"; SRC="$EX/08_icesheet_regional"
    rm -rf "$D/data"; mkdir -p "$D/data"
    cp "$SRC/input/3Dstations.txt" "$D/data/"
    copytree "$SRC/output/stations/JW_stations" "$D/data/stations"
    run_fig "$D"
    ;;
ex09)
    echo "=== ex09 Moon impact — compact surface-field npz (snapshots, record section, seismogram) ==="
    D="$PLOTS_DIR/ex09_moon"; SRC="$EX/09_wave_visualization_Moon/sim"
    rm -rf "$D/data"; mkdir -p "$D/data"
    "$CONDA_AXISEM3D" - "$SRC" "$D/data/moon_surface.npz" <<'PYEOF'
import sys, os, glob, numpy as np, netCDF4
SRC, OUT = sys.argv[1], sys.argv[2]
def readable(c):
    fs = sorted(glob.glob(os.path.join(c, 'axisem3d_synthetics.nc.*')))
    if not fs: return None
    try:
        for f in fs: netCDF4.Dataset(f).close()
        return fs
    except Exception: return None
res = sta = None
for tag in ('1s', '2s', '4s'):   # prefer the 1s run (lunar band; 3000 s record)
    c = os.path.join(SRC, f'output_{tag}', 'stations', 'surface_line')
    if readable(c): res, sta = tag, c; break
if sta is None: sys.exit('  no complete Moon surface_line output found')
lon_of = {}
for L in open(os.path.join(SRC, f'input_{res}', 'moon_line.txt')):
    L = L.strip()
    if not L or L.startswith('#'): continue
    cc = L.split(); lon_of[f'{cc[1]}.{cc[0]}'] = float(cc[3])
time = None; names = []; waves = []
for rf in sorted(glob.glob(os.path.join(sta, 'axisem3d_synthetics.nc.*'))):
    d = netCDF4.Dataset(rf)
    if time is None: time = d.variables['data_time'][:]
    for s in d.variables['list_station'][:]:
        names.append(''.join([c.decode() if isinstance(c, bytes) else str(c)
                              for c in np.ma.filled(s, b'')]).replace('\x00', '').strip())
    waves.append(d.variables['data_wave'][:]); d.close()
wave = np.concatenate(waves, axis=0)
dist = np.array([lon_of.get(n, np.nan) for n in names])
o = np.argsort(dist); dist = dist[o]; wave = wave[o]
ti = np.linspace(0, len(time) - 1, 900).astype(int)
ip = int(np.argmin(np.abs(dist - 90.0)))
np.savez_compressed(OUT, dist=dist.astype('float32'), time=time[ti].astype('float32'),
    uz=wave[:, 2, :][:, ti].astype('float32'), seis_time=time.astype('float32'),
    seis=wave[ip].astype('float32'), seis_dist=float(dist[ip]), res=res)
print('  extracted', res, 'run:', len(dist), 'receivers')
PYEOF
    run_fig_py "$CONDA_AXISEM3D" "$D"
    ;;
ex10)
    echo "=== ex10 handlers demo (50s) — station figures only ==="
    # The full element cross-section (~3.6 GB) is excluded by bundle policy;
    # this bundle ships only the small per-rank station NETCDF group.
    D="$PLOTS_DIR/ex10_handlers"; SRC="$EX/10_handlers_demo/sim_run"
    rm -rf "$D/data"; mkdir -p "$D/data/stations/GSN_Station_Grid"
    cp "$SRC/input/GSN_small.txt" "$SRC/input/inparam.source.yaml" "$D/data/"
    cp "$SRC/output/stations/GSN_Station_Grid/"axisem3d_synthetics.nc.rank* \
       "$SRC/output/stations/GSN_Station_Grid/rank_station.info" \
       "$D/data/stations/GSN_Station_Grid/" 2>/dev/null || echo "    MISSING: ex10 GSN station group"
    run_fig_py "$CONDA_AXISEM3D" "$D"
    ;;
ex11)
    echo "=== ex11 adjoint kernels (1-D PREM + 3-D S362ANI) — axikernels HDF5 ==="
    # Both kernel cases live in ONE bundle (data/kernel_1D + data/kernel_3D); the
    # two source runs sit under examples/11_kernels/{kernel_1D,kernel_3D}.
    D="$PLOTS_DIR/ex11_kernels"; SRC="$EX/11_kernels"
    rm -rf "$D/data"; mkdir -p "$D/data/kernel_1D" "$D/data/kernel_3D"
    # 1-D PREM case
    cp "$SRC/kernel_1D/kernel_output/vp_kernel.h5" "$SRC/kernel_1D/kernel_output/moho_kd.h5" \
        "$D/data/kernel_1D/" 2>/dev/null || echo "    MISSING: kernel_1D HDF5"
    # 3-D S362ANI + Moho-topography case
    cp "$SRC/kernel_3D/kernel_output/vp_kernel.h5" "$SRC/kernel_3D/kernel_output/moho_kd.h5" \
        "$D/data/kernel_3D/" 2>/dev/null || echo "    MISSING: kernel_3D HDF5"
    cp "$SRC/kernel_3D/input_forward/moho_topography.nc" "$D/data/kernel_3D/" \
        2>/dev/null || echo "    MISSING: kernel_3D moho_topography.nc"
    run_fig_py "$CONDA_AXIKERN" "$D"
    ;;
sun)
    echo "=== sun helioseismology — polar map, seismogram, record section ==="
    # ship only the a0 azimuth line plus two mid-distance stations (3-comp seismogram)
    D="$PLOTS_DIR/sun_helioseismology"; SRC="$EX/../../Ben_sun"
    rm -rf "$D/data"; mkdir -p "$D/data/stations"
    cp "$SRC/input/stations_for_plotting.txt" "$SRC/input/inparam.source.yaml" "$D/data/"
    STA="$SRC/output/stations/stationsforplotting"
    cp "$STA/data_time.ascii" "$D/data/stations/" 2>/dev/null || echo "    MISSING: sun data_time.ascii"
    cp "$STA"/a0.*.ascii "$D/data/stations/" 2>/dev/null
    cp "$STA/a30.d60.ascii" "$STA/a60.d60.ascii" "$D/data/stations/" 2>/dev/null || true
    run_fig_py "$CONDA_AXISEM3D" "$D"
    ;;
all)
    for e in ex03 ex07 ex08 ex09 ex10 ex11 sun; do "$0" "$e"; done
    ;;
*)
    echo "Usage: $0 <ex00|ex01|ex02|ex03|ex04|ex05|ex06|ex07|ex08|ex09|ex10|ex11|sun|all>"
    exit 1
    ;;
esac
echo "=== package $1 done ==="
