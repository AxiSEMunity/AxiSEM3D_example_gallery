# Example 3: SEG/EAGE salt body, local Cartesian mesh (`03_salt_body_SEG_local`)

SEG/EAGE C3 salt model (Aminzadeh et al., 1997) in a Cartesian box of 5 km radius with a 200-m water
layer and an undulating sea floor; moment tensor source on the axis at 1 km depth; 10 s. Shown in
Fig. 3 and Table 2 (row 3).

| | |
|---|---|
| Mesh period | 0.2 s (5 Hz) |
| Mesh | `input3D/local_mesh__SEG_salt__5Hz.e`, 2,000 elements |
| 3-D models | `SEG_C3_data/SEG_C3_SOLID.nc` (volumetric) and `SEG_C3_data/SEG_C3_UNDULATION.nc` (sea floor), both packed in `input3D/SEG_C3_data.tar.bz2` |
| Azimuthal expansion | `type_Nr: CONSTANT`, 50 |
| Source time function | Gaussian, half duration 0.2 s |
| Output | element groups `orthogonal_azimuthal_slices` and `Fourier_coefficients_ocean_floor` |
| MPI ranks used | 4 |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEMCartesian --basic.model SEG_C3.bm --basic.period .2 --cartesian2Daxisem.x 5. --cartesian2Daxisem.min_z 6367.0 --attenuation.frequencies 0.01 10. --output_filename local_mesh__SEG_salt__5Hz.e
```

## Run

```
mkdir -p input3D/SEG_C3_data && tar -xjf input3D/SEG_C3_data.tar.bz2 -C input3D/SEG_C3_data
mpirun -np 4 axisem3d --input input3D --output output3D
```

The archive unpacks to `SEG_C3_SOLID.nc` (1.05 GB) and `SEG_C3_UNDULATION.nc` (1.8 MB).

## Files

- `input3D/`: the five `inparam.*.yaml` files, mesh, `SEG_C3_data.tar.bz2` (43 MB).
- `SEG_C3.bm`: 1-D background model of the mesh.
- `make_figures.py`: figures from `output/`.
- `paper_figure/make_figures.py`, `pub_style.py`: script of the panels of Fig. 3. It reads `data/ocean_floor`
  and `data/azimuthal_slices` (copies of the two element-output folders) and `data/salt_model.npz` (two slices
  of `SEG_C3_SOLID.nc`); see case `ex03` in `../paper_figure_tools/package.sh`.

No file is larger than 50 MB. The inputs are those of `examples/03_salt_body_SEG_local/input3D` of tag `v2.1.0`,
unchanged.
