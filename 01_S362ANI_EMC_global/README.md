# Example 2a: global wavefield with 3-D mantle structure (`01_S362ANI_EMC_global`)

Anisotropic PREM mesh with the mantle model S362ANI (Kustowski et al., 2008) added as a volumetric
3-D model; same source and stations as Example 1; 1800 s. Shown in Fig. 2b and Table 2 (row 2a).

| | |
|---|---|
| Mesh period | 5 s |
| Mesh | `global_prem_ani_5s.e`, 478,222 elements (the mesh of Example 1) |
| 3-D model | `S362ANI_percent.nc` (IRIS EMC file; `dvs` in percent, applied to VS with factor 0.01 and to VP with factor 0.005, relative to the 1-D model) |
| Azimuthal expansion | `type_Nr: CONSTANT`, 50 |
| Source time function | Gaussian, half duration 2.5 s |
| Stations | GSN (129); USArray TA (1,834; 1 s sampling); 0.1-degree surface grid (`US_grid_0p1deg.txt`, 208,681 points; 1 s sampling) |
| Courant number | 0.5 |
| MPI ranks used | 384 |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani.bm --basic.period 5 --output_filename global_prem_ani_5s.e
```

## Run

```
mpirun -np 384 axisem3d --input input_5s --output output_5s
```

## Files

- `input_5s/`: the five `inparam.*.yaml` files, station files, `S362ANI_percent.nc`, mesh.
- `input/S362ANI_percent.nc`: the same model file, for the 50-s template of the solver repository.
- `make_figures.py`: figures from `output/`. The script of Fig. 2 is in `../00_global_1D/paper_figure/`.

No file is larger than 50 MB.

## Changes relative to the files of the run

`exodus_mesh:` in `inparam.model.yaml` gave an absolute path on the cluster; it now gives the file name
of the mesh, which is placed in the input folder. Nothing else was changed.
