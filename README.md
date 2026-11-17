# AxiSEM3D example gallery: inputs of the paper figures

This gallery holds, for each of the eleven examples of Fernando et al., *AxiSEM3D: efficient 3D wave
propagation across scales and worlds* (Geophysical Journal International), the complete solver inputs
of the runs that produced the figures and the entries of Table 2, together with the plotting scripts.

The low-resolution templates of the same examples are distributed with the solver
(<https://github.com/AxiSEMunity/AxiSEM3D>, folder `examples/`). The folders here have the same names.

## Solver and tools

- **Solver:** AxiSEM3D v2.1.0. The runs used commit `ca260c7` of the solver repository; its solver
  source (`src/`) is identical to tag `v2.1.0`.
- **Meshes:** `salvus_mesh_lite` (<https://github.com/AxiSEMunity/salvus_mesh_lite>). Every Exodus mesh
  stores the command that generated it; the README of each example repeats it.
- **Kernels (Example 7):** `axikernels` (<https://github.com/Adrian-Mag/AxiSEM3D_Kernels>), commit `a1c853e`.

## Examples

| Paper | Folder | Shown in | Inputs of the run | Period |
|---|---|---|---|---|
| 1 | `00_global_1D` | Fig. 2a; Table 2 | `input_5s/`; `input_2s/` | 5 s; 2 s |
| 2a | `01_S362ANI_EMC_global` | Fig. 2b; Table 2 | `input_5s/` | 5 s |
| 2b | `02_3d_crust_S362ANI_regional` | Fig. 2c; Table 2 | `simu_with_3d_crust_5s/input/` | 5 s |
| 3 | `03_salt_body_SEG_local` | Fig. 3; Table 2 | `input3D/` | 0.2 s |
| 4 | `04_simple_3d_shapes` | Fig. 4a, 4b; Table 2 | `example_input_cartesian/`, `example_single_plume/`; `example_release_paper/input/` | 0.2 s, 14 s; 0.08 s |
| 5 | `05_anisotropy_global` | Fig. 5; Table 2 | `2012-07-03_paper_example_50s/sim_US32_olivineE_5s/input/`, `.../sim_US32_isotropic_5s/input/` | 5 s |
| 6 | `06_finite_source_local` | Fig. 6; Table 2 | `input/` (model files: Zenodo record of the model, see README) | 2 s |
| 7 | `07_kernels` | Fig. 7; Table 2 | `kernel_1D/`, `kernel_3D/` (forward and adjoint), `fig7d_test/` | 50 s |
| 8 | `08_atmosphere_Mars_global` | Fig. 8a; Table 2 | `axisem3d_mars_atm/input/` | 20 s |
| 9 | `09_icesheet_regional` | Fig. 8b; Table 2 | `input/` | 2 s |
| 10 | `10_helioseismology_sun` | Fig. 8c; Table 2 | `input/` | 500 s |
| 11 | `11_wave_visualization_Moon` | Fig. 9; Table 2 | `sim/input_scatter/` | 2 s |

Examples 7 and 10 were added to the solver repository after tag `v2.1.0`; their templates are on its
`main` branch.

## What an example folder contains

- `README.md`: what the run is, period, mesher command, run command, figure, large files.
- The input folder(s) of the run: the five `inparam.*.yaml` files, station lists, the mesh, and the
  3-D model files. An input folder is complete: `axisem3d --input <folder>` needs nothing else.
- Plotting scripts in the example folder (`make_figures.py` and others). They read the `output/`
  folder of a run.
- `paper_figure/`: the script(s) as used for the figure in the paper. Most of them read compact
  extracts of the solver output from a folder `data/` next to the script. `paper_figure_tools/package.sh`
  records how each extract was made from the output (its folder names are those of the original run
  tree, listed in the README of each example).

No simulation output is stored here.

## Large files

Files larger than 50 MB are not in the Git repository. They are in the Zenodo archive of this gallery
(the DOI will be added here), under the same relative paths. `MANIFEST.csv` lists every file with its size,
SHA-256 checksum and location.

## Running an example

```
mpirun -np N axisem3d --input <input folder> --output <output folder>
```

The README of each example gives the folder names and the number of MPI ranks that were used.
