# Example 7: sensitivity kernels (`07_kernels`)

Traveltime sensitivity kernel of the direct P wave to v_p for one source-receiver pair 40 degrees apart (source at
0N 0E, 20 km depth; receiver at 0N 40E), from a forward and an adjoint (backward) AxiSEM3D run, and its test against
direct simulation. Shown in Fig. 7 and Table 2 (rows 7): `kernel_1D/` is panel (a), `kernel_3D/` panel (b),
`fig7d_test/` panel (d).

| | `kernel_1D/` | `kernel_3D/` |
|---|---|---|
| Mesh period | 50 s | 50 s |
| Mesh | `input_forward/prem_iso_elastic_50s.e`, 6,400 elements | the same mesh |
| 3-D models | none | `S362ANI_radius.nc` (mantle; VP and VS perturbations multiplied by 5) and `moho_topography.nc` (Gaussian Moho undulation, peak 5 km) |
| Azimuthal expansion | `type_Nr: CONSTANT`, 5 | `type_Nr: CONSTANT`, 50 |
| Forward source | explosion, Gaussian half duration 25 s, 600 s | the same |
| Adjoint source | vertical force at the receiver (`backward_input_as_run/STF/`), 663 s | the same |
| Element output | all GLL points of every element (`GLL_points_one_edge: FULL`) | the same |
| MPI ranks used | 8 | 8 |

## Measurement and kernel

Cross-correlation traveltime of P on the vertical displacement at the receiver, window 425-475 s with 15-s cosine
tapers. The relative v_p kernel is `K = 2 rho vp^2 K_lambda / tau` with `K_lambda = int div s . div s~ dt`
(`dT = int K dln(vp) dV`). The divergence is a scalar, so the forward and adjoint fields need no rotation between the
two run frames. The scripts are in `paper_figure/`; usage is at the top of each:

1. `make_adjoint_source.py`: adjoint source from the forward receiver record (writes `STF/` of the backward input).
2. `eval_kvp.py`: K_lambda and the kernel on a slice and on the test perturbation (reads the element output of both runs
   with `axk_full.py`, a reader for all GLL points built on `axikernels`).
3. `test_kernels.py`: Fig. 7d, kernel prediction against the traveltime change measured in the direct runs.
4. `plot_fig7.py`: the figure.

The scripts take the runs from one folder (`KERNEL_WORK`): `redo_full/{fwd,bwd}{1d,3d}` = the runs of
`kernel_*/input_forward` and `kernel_*/backward_input_as_run`; `redo/` holds `adjoint_source_<dim>.npz` and the forward
receiver record (`redo/fwd<dim>`, which may be a link to `redo_full/fwd<dim>`: the receiver records are the same);
`runs/ref<dim>` and `runs/<dim>_B_*` = the runs of `fig7d_test/`.

## Fig. 7d: direct simulations

`fig7d_test/<dim>_B_{p050,p100,m050,m100}/input/`: the forward setup (receiver output only) with a Gaussian v_p
perturbation (`blob_vp.nc`, 1 % at the center) scaled by +0.5, +1, -0.5 and -1; `fig7d_test/ref<dim>/input/`: the same
without the perturbation. 8 MPI ranks each.

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM -o prem_iso_elastic_50s.e --overwrite_file --basic.model prem_iso_elastic.bm --basic.period 50
```

## Template

`kernel_*/compute_kernels*.py`, `run.sh` and `readme.txt` are the template workflow of this example with `axikernels`
(<https://github.com/Adrian-Mag/AxiSEM3D_Kernels>, commit `a1c853e`). Fig. 7 was made with the scripts in
`paper_figure/` instead.

## Files

- `kernel_*/input_forward/`, `kernel_*/backward_input_as_run/`: the five `inparam.*.yaml` files, station files,
  meshes, `prem_iso_elastic.bm`, (3-D) the two model files, and (backward) the adjoint source `STF/`.
- `kernel_3D/convert_s362ani_to_radius.py`, `create_moho_topography.py`: scripts that built the two model files.
- `fig7d_test/`: inputs of the ten direct runs of Fig. 7d.

No file is larger than 50 MB. The inputs were not changed relative to the files of the runs; `fig7d_test/ref*/input/`
were assembled from the files of the test runs and the forward runs (identical checksums to the reference runs).
