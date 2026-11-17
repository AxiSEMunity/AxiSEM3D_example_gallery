# Example 6: kinematic rupture in the San Francisco Bay Area (`06_finite_source_local`)

Mw 7.1 kinematic rupture on the Hayward Fault in the 3-D USGS San Francisco Bay Area model, on a local Cartesian
mesh with surface topography and no velocity floor (minimum VS 80 m/s). Shown in Fig. 6 and Table 2 (row 6).

| | |
|---|---|
| Mesh period | 2 s |
| Mesh | `AxiSEMCartesian_sfba_2018_2s.e`, 33,900 elements; a cylinder of 50 km radius whose axis passes through the model center (37.85517N, 122.25684W) |
| 1-D background | `sfba_2018.bm` |
| 3-D model | four velocity layers, surface topography and Moho (see "Large files") |
| Azimuthal expansion | `SFBA_finite_rupture_Nr0.nc` (maximum N_r 9,720; sum over the mesh nodes 13.1 million) |
| Source | 2,647 sub-sources in `inparam.source.yaml` (moment tensors in dyne cm, `unit: 1e-07`; total moment 5.78e19 N m, Mw 7.1; Gaussian slip-rate functions, half duration 2 s, with rupture-time shifts) |
| Attenuation | `FULL` |
| Time step | Courant number 0.5 |
| Record length | 25 s |
| Receivers | `STATIONS_OUTPUT.txt`: 101,765 surface receivers 250 m apart within 45 km of the axis |

**Stability.** With the coarse-grained attenuation (`attenuation: CG4`) this run becomes unstable after about 24 s
near the mesh axis, where the 3-D model has a 25-m surface layer with VS 120 m/s and Q_mu 13. With `FULL` it is
stable over the 25-s record. Fig. 6 shows t = 16 s and the peak ground velocity up to 16 s.

## Source list

The sub-source list shipped with AxiSEM3D v2.1.0 (`examples/06_finite_source_local/input/inparam.source.yaml` of the
tag, 36,246 sub-sources) holds its moment tensors in one fixed (vertical, north, east) frame, while the solver reads
the tensor of a source given by distance and azimuth in the local (z, s, phi) frame at that azimuth.
`source_list/build_rupture_2s.py` rotates every tensor into the solver's frame and merges the sub-sources on a 400-m
grid along strike and depth (tensors summed; position and rupture time moment-weighted), which leaves the total moment
tensor unchanged. `source_list/fix_source_locations.py` then moves 14 sub-sources toward the center of their element
(median 47 m), so that the solver can locate them. The header of `inparam.source.yaml` repeats this.

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model sfba_2018.bm --basic.period 2.0 --attenuation.frequencies 0.01 10.0 --cartesian2Daxisem.x 50.0 --output_file AxiSEMCartesian_sfba_2018_2s.e
```

## Run

```
mpirun -np 512 axisem3d --input input --output output
```

## Fig. 6

```
HDF5_USE_FILE_LOCKING=FALSE python paper_figure/fig6_rupture.py . ex06_wavefield_snapshots --tsnap 16 --pgvtmax 16
```

Shoreline: `paper_figure/coast_sfba_gshhg_full_L1.txt` (GSHHG full resolution, level 1).

## Large files: existing Zenodo record

The five large model files are the members of `SFBA.zip` in <https://doi.org/10.5281/zenodo.20017355>
(Haindl and Muir), under other names. Place them in `input/` as:

| Name in `input/` | Name in `SFBA.zip` | Size |
|---|---|---|
| `SFBA_layer1_min500.nc` | `SFBA_layer1_2018.nc` | 761 MB |
| `SFBA_layer2.nc` | `SFBA_layer2_2018.nc` | 452 MB |
| `SFBA_layer3.nc` | `SFBA_layer3_2018.nc` | 67 MB |
| `SFBA_layer4.nc` | `SFBA_layer4_2018.nc` | 45 MB |
| `topo_smooth_100m_upsampled.nc` | `topo_smooth_100m_upsampled_2018.nc` | 401 MB |

`SFBA_layer1_min500.nc` holds no velocity floor despite its name (it is identical to `SFBA_layer1_2018.nc`).
`moho_smooth_24km.nc` (in `input/`) is `moho_smooth_24km_2018.nc` of the same archive.

## Changes relative to the files of the run

`inparam.source.yaml`: `record_length` 25 s (the run behind Fig. 6 used 40 s; the first 25 s are the same), and the
first comment line no longer names the manuscript. Nothing else was changed.
