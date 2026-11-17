# Example 4: simple 3-D shapes (`04_simple_3d_shapes`)

Three runs. Each folder is a complete input folder (the `inparam.*.yaml` files are in the folder itself).

| | `example_input_cartesian/` | `example_single_plume/` | `example_release_paper/input/` |
|---|---|---|---|
| Shown in | Fig. 4a | Fig. 4b | Table 2, row 4 (not illustrated) |
| Structure | one sphere in a homogeneous Cartesian box | plume (head and conduit) on the axis of a homogeneous sphere of radius 4000 km | homogeneous Cartesian box with 1,000 random ellipsoids in one quadrant and a regular lattice of spheres in another |
| Mesh period | 0.2 s | 14 s | 0.08 s |
| Mesh | `homogenous_cartesian.e`, 42,670 elements | `homogenous_global.e`, 20,178 elements | `paper_example.e`, 281,250 elements |
| 3-D model | `example_sphere_cartesian.nc` (VP -20 %, VS -100 %, density -20 % inside the sphere) | `plume_1.nc` (VP, VS, density -10 %) | `paper_example.nc` (VP, VS, density; perturbations down to -47 %) |
| Azimuthal expansion | `CONSTANT`, 360, with wavefield scanning | `CONSTANT`, 5 | `CONSTANT`, 360 |
| Source | moment tensor on the axis, 8 km depth, Gaussian half duration 0.2 s (the mesh period) | as 4a, half duration 14 s (the mesh period) | moment tensor on the axis, 7.5 km depth, Gaussian half duration 0.5 s |
| Courant number | 0.4 | 0.5 | 0.4 |
| Record length | 15 s | 1200 s | 50 s |
| MPI ranks used | 64 | 4 | 144 |

The model files hold fractional perturbations relative to the 1-D model (`reference_kind: REF1D`).

## Meshes

```
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model homogenous_cart.bm --basic.period .2 --cartesian2Daxisem.x 25. --cartesian2Daxisem.min_z 0.0 --attenuation.frequencies 0.01 10. --output_filename homogenous_cartesian.e
python -m salvus_mesh_lite.interface AxiSEM --basic.model homogenous_global.bm --basic.period 14 --output_filename homogenous_global.e
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model homogenous.bm --basic.period .08 --cartesian2Daxisem.x 30. --cartesian2Daxisem.min_z 0.0 --attenuation.frequencies 0.01 14. --output_filename paper_example.e
```

These are the commands stored in the three mesh files. `example_release_paper/input/gen_mesh.sh` gives
`--attenuation.frequencies 0.01 16.`; the mesh of the run was made with `0.01 14.`.

## 3-D models

`example_input_cartesian/gen_blob.py`, `example_single_plume/gen_plume.py` and
`example_release_paper/generate_3D_model.py` build the model files with the classes in `src/`.
`generate_3D_model.py` reads the 1,000 random ellipsoids from `example_release_paper/randomblobs` and adds the regular lattice.
The model files of the runs are included, so the scripts need not be run.

## Run

```
mpirun -np 64  axisem3d --input example_input_cartesian     --output example_input_cartesian/output
mpirun -np 4   axisem3d --input example_single_plume        --output example_single_plume/output
mpirun -np 144 axisem3d --input example_release_paper/input --output example_release_paper/output
```

## Figures

- `make_figures.py`, `make_blob_plume_figure.py`: figures from the outputs.
- `paper_figure/make_blob_plume_figure.py`, `pub_style.py`: script of Fig. 4. It reads compact files in
  `data/` (model slices and wavefield slices at azimuth 0), written from the run outputs by
  `paper_figure/package_ex4a.py` and `package_ex4b.py` (usage at the top of each).
- `example_release_paper/reproduce_paper_figure.py`: record section of the 1,000-ellipsoid run.

## Large files (Zenodo archive, not in the Git repository)

- `example_input_cartesian/example_sphere_cartesian.nc` (73 MB)
- `example_single_plume/plume_1.nc` (75 MB)
- `example_release_paper/input/paper_example.nc` (5.18 GB)

The inputs were not changed relative to the files of the runs.
