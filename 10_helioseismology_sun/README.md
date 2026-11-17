# Example 10: the Sun (`10_helioseismology_sun`)

Acoustic waves in a 1-D, fully fluid solar model, excited by a pressure source; 12,000 s. Shown in Fig. 8c
and Table 2 (row 10).

| | |
|---|---|
| Mesh period | 500 s, with 4 elements per wavelength |
| Mesh | `input/solar_master_epw4.e`, 1,080,586 elements, all fluid |
| 3-D models | none (the model listed in `inparam.model.yaml` is deactivated) |
| Azimuthal expansion | `type_Nr: CONSTANT`, 7 |
| Source | `FLUID_PRESSURE` at latitude 89.9, longitude 45, depth 5.0e7 m; Gaussian, half duration 382 s (the mesh period, 381.6 s, so that the source holds no content the mesh cannot resolve) |
| Time scheme | `SYMPLECTIC`, Courant number 0.4 |
| Output | 57 stations (`stations_for_plotting.txt`, 0.05 s sampling); element group `sun_cross_section` on the plane of azimuth 0 (50 s sampling) |
| MPI ranks used | 896 |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model solar_master.bm --basic.period 500 --advanced.elements_per_wavelength 4 --output_filename solar_master_epw4.e
```

## Run

```
mpirun -np 896 axisem3d --input input --output output
```

## Files

- `input/`: the five `inparam.*.yaml` files, `stations_for_plotting.txt`, `solar_master.bm`, mesh.
- `input/solar_new_test_500.e`: smaller mesh (272,160 elements) for a quick test; not used for the paper.
- `make_figures.py`, `make_cross_section.py`: figures from `output/`.
- Fig. 8 as printed: `../paper_figure_tools/dump_sun_v6.py` (snapshots from `output/elements/sun_cross_section/`) and
  `../paper_figure_tools/fig8_composite.py`, mode `v6` (records low-passed at the mesh period).
- `paper_figure/make_cross_section.py`, `pub_style.py`: script of Fig. 8c as submitted. It reads `output/elements/sun_cross_section/`
  and `output/stations/stationsforplotting/`. For the paper it was run through
  `../paper_figure_tools/bigfont_runner.py` with `FONT_SCALE=1.5`.

No file is larger than 50 MB.

## Changes relative to the files of the run

None in the inputs. In `paper_figure/make_cross_section.py` the absolute cluster path `RUN` was replaced by
the path relative to this folder.
