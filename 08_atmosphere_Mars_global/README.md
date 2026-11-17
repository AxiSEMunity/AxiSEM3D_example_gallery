# Example 8: Mars with an atmosphere (`08_atmosphere_Mars_global`)

Pressure source in a 60-km fluid atmosphere on top of the TAYAK interior model; 1800 s. Shown in Fig. 8a
and Table 2 (row 8).

| | |
|---|---|
| Mesh period | 20 s |
| Mesh | `axisem3d_mars_atm/input/tayak_60km_20s.e`, 204,456 elements |
| 3-D models | none |
| Azimuthal expansion | `type_Nr: CONSTANT`, 1 |
| Source | `FLUID_PRESSURE` at radius 3,394.0 km (the solid surface is at 3,389.5 km), half duration 0 (impulse) |
| Absorbing boundary | top of the atmosphere (Clayton-Enquist and sponge, relative span 0.01) |
| Output | 10 surface stations (`mars.txt`, 0 to 9 degrees); element group `mars_xs_fluid`: the atmosphere on the plane of azimuth 0, 0 to 7 degrees, 0.5 s sampling (so that it can be low-passed at the 20-s mesh period without aliasing) |
| MPI ranks used | 128 |

## Mesh

```
cd axisem3d_mars_atm && bash gen_mesh.sh 20
```

`gen_mesh.sh` rewrites the atmosphere of `tayak_60km_atm.bm` on a denser radius grid (`tayak_60km_atm_meshing.bm`),
runs `python -m salvus_mesh_lite.interface AxiSEM --basic.model tayak_60km_atm_meshing.bm --basic.period 20
--output_filename input/tayak_60km_20s.e`, and checks that the fluid density of the mesh is positive.

## Run

```
mpirun -np 128 axisem3d --input axisem3d_mars_atm/input --output axisem3d_mars_atm/output_xs
```

## Files

- `axisem3d_mars_atm/input/`: the five `inparam.*.yaml` files, `mars.txt`, mesh.
- `axisem3d_mars_atm/`: `gen_mesh.sh`, the two `.bm` files, `prepare_positive_atmosphere_bm.py`,
  `validate_exodus_fluid_rho.py`, `README.txt` (notes on atmosphere models and the absorbing boundary).
- `make_figures.py`: figures from the station output.
- Fig. 8 as printed (rows a to c, low-passed at the mesh periods): `../paper_figure_tools/fig8_composite.py`, mode `v6`.
- `paper_figure/make_atm_section.py`, `pub_style.py`: script of Fig. 8a as submitted. It reads the element group
  `mars_xs_fluid` of `axisem3d_mars_atm/output_xs/` and, for the record section, `data/mars.txt` and
  `data/stations/` (copy of the station output; case `ex07` in `../paper_figure_tools/package.sh`, where this
  example is the folder `07_atmosphere_Mars_global`). For the paper it was run through
  `../paper_figure_tools/bigfont_runner.py` with `FONT_SCALE=1.5`.

No file is larger than 50 MB.

## Changes relative to the files of the run

`inparam.output.yaml`: the comment above the element group was shortened. In `paper_figure/make_atm_section.py` the absolute cluster path `RUN` was replaced by
the path relative to this folder.
