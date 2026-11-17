# Example 1: global wavefield in a 1-D model (`00_global_1D`)

Anisotropic PREM, 2011 Mw 5.8 Virginia earthquake (moment tensor source at 12 km depth), 1800 s.

| | `input_5s/` | `input_2s/` |
|---|---|---|
| Shown in | Fig. 2a | Table 2, row 1 |
| Mesh period | 5 s | 2 s |
| Mesh | `global_prem_ani_5s.e`, 478,222 elements | `global_prem_ani_2s.e`, 2,911,386 elements |
| Azimuthal expansion | `type_Nr: CONSTANT`, 5 | `type_Nr: CONSTANT`, 5 |
| Source time function | Gaussian, half duration 2.5 s | Gaussian, half duration 0 (impulse) |
| Stations | GSN (129); USArray TA (1,834; 1 s sampling); 0.1-degree surface grid over the United States (`US_grid_0p1deg.txt`, 208,681 points; 1 s sampling) | GSN (129), USArray TA (1,834; 5 s sampling) |
| Courant number | 0.5 | 0.5 |
| MPI ranks used | 128 | 216 |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani.bm --basic.period 5 --output_filename global_prem_ani_5s.e
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani.bm --basic.period 2 --output_filename global_prem_ani_2s.e
```

`prem_ani.bm` is the anisotropic PREM file that ships with `salvus_mesh_lite`.

## Run

```
mpirun -np 128 axisem3d --input input_5s --output output_5s
mpirun -np 216 axisem3d --input input_2s --output output_2s
```

## Fig. 2

`paper_figure/make_fig2.py` draws all three rows of Fig. 2 (Examples 1, 2a and 2b) from the station group
`US_grid_0p1deg` of the three 5-s runs: vertical displacement maps, zero-phase low-pass filtered at 15 s.
Shorelines and borders are in `paper_figure/map_lines/` (GSHHG level 1 and GMT borders, dumped with
`dump_lines.sh`). Usage is at the top of the script.

## Files

- `input_5s/`, `input_2s/`: `inparam.{model,nr,source,output,advanced}.yaml`, station files, mesh.
- `make_figures.py`: seismogram, station map and USArray time-lapse figures from `output/`.
- `paper_figure/`: the script of Fig. 2 and its map lines.

**Large file:** `input_2s/global_prem_ani_2s.e` (172 MB) is in the Zenodo archive, not in the Git repository.

## Changes relative to the files of the run

`exodus_mesh:` in `inparam.model.yaml` gave an absolute path on the cluster; it now gives the file name
of the mesh, which is placed in the input folder. Nothing else was changed.
