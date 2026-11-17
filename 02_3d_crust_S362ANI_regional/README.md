# Example 2b: 3-D crust and mantle on a regional mesh (`02_3d_crust_S362ANI_regional`)

Anisotropic PREM mesh limited to 90 degrees from the source, with S362ANI in the mantle and CRUST1.0
(Laske et al., 2013) in the crust; same source as Example 1; 1800 s. Shown in Fig. 2c and Table 2 (row 2b).

| | |
|---|---|
| Mesh period | 5 s |
| Mesh | `regional_prem_ani_5s.e`, 194,994 elements |
| 3-D models | `S362ANI_percent.nc` (mantle); `crust1_volumetric.nc` (VP, VS, density of the crust, absolute values); `crust1_geometric.nc` (surface and Moho undulation, each applied with factor 0.5); `crust1_oceanload.nc` (ocean load) |
| Azimuthal expansion | `type_Nr: CONSTANT`, 100 |
| Source time function | Gaussian, half duration 2.5 s |
| Stations | USArray TA (1,834; 1 s sampling); 0.1-degree surface grid (`US_grid_0p1deg.txt`, 208,681 points; 1 s sampling) |
| Courant number | 0.6 |
| MPI ranks used | 256 |

CRUST1.0 enters as three external structured-grid files (classes `StructuredGridV3D`, `StructuredGridG3D`,
`StructuredGridO3D`), as in the template of tag `v2.1.0`.

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_ani.bm --basic.period 5 --chunk2D.max_colatitude 90 --output_filename regional_prem_ani_5s.e
```

## Run

```
mpirun -np 256 axisem3d --input simu_with_3d_crust_5s/input --output simu_with_3d_crust_5s/output
```

## Files

- `simu_with_3d_crust_5s/input/`: the five `inparam.*.yaml` files, station files, the four model files, mesh.
- `gen_crust1.ipynb`: notebook that built the three `crust1_*.nc` files from the CRUST1.0 distribution
  (the raw CRUST1.0 files are in `examples/02_3d_crust_S362ANI_regional/crust1.0_raw_data/` of tag `v2.1.0`).
- `input_share/S362ANI_percent.nc`, `input_with_3d_crust/crust1_volumetric.nc`: the same model files, for
  the 50-s template of the solver repository.
- `make_figures.py`: figures from the output. The script of Fig. 2 is in
  `../00_global_1D/paper_figure/`.

No file is larger than 50 MB.

## Changes relative to the files of the run

`exodus_mesh:` in `inparam.model.yaml` gave an absolute path on the cluster; it now gives the file name
of the mesh, which is placed in the input folder. Nothing else was changed.
