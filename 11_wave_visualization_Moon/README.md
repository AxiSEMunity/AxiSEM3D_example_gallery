# Example 11: impact on the Moon with 3-D scattering (`11_wave_visualization_Moon`)

Isotropic source at 1 km depth in the lunar model ISSI_MOON_M1 (Garcia et al., 2019) with a 3-D stochastic
scattering model in the uppermost 55 km; 1000 s. Shown in Fig. 9 and Table 2 (row 11). The inputs of this
run are in `sim/input_scatter/`.

| | |
|---|---|
| Mesh period | 2 s |
| Mesh | `sim/input_scatter/moon_issi_m1_2s.e`, 1,001,664 elements |
| 3-D model | `sim/input_scatter/moon_scatter.nc`: fractional perturbations of VP, VS and density relative to the 1-D model, on a grid of 1 degree by 1 degree by 5 km |
| Azimuthal expansion | `type_Nr: CONSTANT`, 256 |
| Source | moment tensor (explosion) at 0N 0E, Gaussian half duration 2.5 s |
| Output | `moon_grid.txt`: 7,080 surface receivers on a 3-degree grid (vertical displacement, 3 s sampling); `moon_line.txt`: 1,799 receivers along the equator |
| MPI ranks used | 1,024 |

The scattering model: the VP perturbation has an rms of 9 % at the surface and 6 % between 15 and 40 km depth
and is tapered to zero at 60 km. The VS perturbation is 1.25 times and the density perturbation 0.625 times the
VP perturbation at every grid point.

## Mesh

```
cd sim && python -m salvus_mesh_lite.interface AxiSEM --basic.model ISSI_MOON_M1_atten.bm --basic.period 2 --output_filename input_scatter/moon_issi_m1_2s.e
```

## Run

```
mpirun -np 1024 axisem3d --input sim/input_scatter --output sim/output_scatter
```

## Files

- `sim/input_scatter/`: the five `inparam.*.yaml` files, `moon_grid.txt`, `moon_line.txt`, mesh, `moon_scatter.nc`.
- `sim/ISSI_MOON_M1_atten.bm`: 1-D background model of the meshes.
- `sim/input/moon_issi_m1_{1s,2s,4s}.e`: meshes for the 1-D, axisymmetric template of the solver repository;
  not used for the figure.
- `paper_figure/make_3d_snapshots_scatter_hybrid_tight.py`, `pub_style.py`: script of Fig. 9. It reads
  `sim/input_scatter/moon_grid.txt`, the station group `surface_grid` of `sim/output_scatter/`, and
  `assets/moon_8k.jpg`.
- `make_3d_snapshots_scatter.py` and the other `make_*.py` scripts, `FIGURES.md`: further figures and movies.
- `assets/`: lunar surface texture used as background.

**Large file:** `sim/input/moon_issi_m1_1s.e` (158 MB) is in the Zenodo archive, not in the Git repository.

## Changes relative to the files of the run

None in the inputs. In `paper_figure/make_3d_snapshots_scatter_hybrid_tight.py` the absolute cluster path `RUN`
and the two texture paths were replaced by paths relative to this folder.
