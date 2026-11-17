# Example 5: full anisotropy, SKS splitting (`05_anisotropy_global`)

Deep earthquake of 3 July 2012 beneath Cook Strait, New Zealand, recorded at 180 stations in the
northwestern United States; 2000 s. Two runs on the same isotropic PREM mesh: one with an olivine
elastic tensor in the upper mantle oriented after the US32 model, one isotropic reference.
Shown in Fig. 5 and Table 2 (row 5).

All paths below are inside `2012-07-03_paper_example_50s/`.

| | `sim_US32_olivineE_5s/input/` | `sim_US32_isotropic_5s/input/` |
|---|---|---|
| Role in Fig. 5 | anisotropic synthetics (blue) and splitting measurements | isotropic reference (black) |
| Mesh period | 5 s | 5 s |
| Mesh | `global_prem_iso_5s.e`, 475,710 elements | the same mesh |
| 3-D model | `jons_azi_410_3.nc`: 21 elastic constants (GPa), absolute values | none (the model entry is deactivated) |
| Azimuthal expansion | `type_Nr: ANALYTICAL`, 100 at all depths | `type_Nr: CONSTANT`, 5 |
| Source | moment tensor, depth 229.8 km, half duration 0 (impulse) | the same |
| Stations | `eq3.txt` (180) | the same |
| MPI ranks used | 216 | 288 |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model prem_iso.bm --basic.period 5 --output_filename global_prem_iso_5s.e
```

`prem_iso.bm` is the isotropic PREM file that ships with `salvus_mesh_lite`.

## Run

```
mpirun -np 216 axisem3d --input sim_US32_olivineE_5s/input  --output sim_US32_olivineE_5s/output
mpirun -np 288 axisem3d --input sim_US32_isotropic_5s/input --output sim_US32_isotropic_5s/output
```

## Files

- The two input folders: the five `inparam.*.yaml` files, `eq3.txt`, mesh, and (anisotropic run) `jons_azi_410_3.nc`.
- `processing/US32_olivine_410.m`: MATLAB script that built the elastic-tensor model from US32.
- `sim_US32_olivineE/input/jons_azi_410_3.nc`: the same model file, for the 50-s template of the solver repository.
- `../make_figures.py`, `../make_composite_figure.py`, `../make_splitting_figure.py`: figures from the outputs.
- `../paper_figure/make_fig5.py`: script of Fig. 5 (both panels; records low-passed at the 5-s mesh period).
  It reads the station output of the two runs (`sim_US32_olivineE_5s/output/`, `sim_US32_isotropic_5s/output/`)
  and the files in `../paper_figure/data/`. Usage at the top of the script.
- `../paper_figure/make_composite_figure.py`, `make_splitting_figure.py`, `pub_style.py`: scripts of Fig. 5 as submitted.
  They read the station output of the two runs from `data/sim_US32_olivineE_5s/` and `data/sim_US32_isotropic_5s/`
  (case `ex05` in `../paper_figure_tools/package.sh`) and the files in `../paper_figure/data/`:
  - `us32_stations.txt`: copy of `eq3.txt`;
  - `splitting/20120703_phi_dt_SKS_vel.txt`: fast directions and delay times measured on the anisotropic
    synthetics with SplitRacer (not a solver input; measured outside this repository);
  - `splitting/Long_splitting.txt`: observations of Long et al. (2009).

No file is larger than 50 MB.

## Changes relative to the files of the runs

`exodus_mesh:` in `inparam.model.yaml` gave an absolute path on the cluster; it now gives the file name
of the mesh, which is placed in each input folder. Nothing else was changed.
