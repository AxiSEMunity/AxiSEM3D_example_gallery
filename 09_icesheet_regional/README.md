# Example 9: Enceladus ice shell and ocean, regional mesh (`09_icesheet_regional`)

Regional mesh of Enceladus (60 degrees from the source, down to 150 km radius) with ice shell, ocean and
silicate interior; moment tensor source at 2 km depth. Shown in Fig. 8b and Table 2 (row 9).

| | |
|---|---|
| Mesh period | 2 s |
| Mesh | `input/enceladus_2s_chunk.e`, 9,856 elements |
| 3-D models | `volumetric_ice_shell_properties.nc` (VP and VS of the ice shell, absolute values; Dapré et al.); `Cadek_basal_topography.nc` (geometric model, variable `iob`; Čadek et al., 2016) |
| Azimuthal expansion | `type_Nr: CONSTANT`, 5 |
| Source time function | half duration 0 (impulse) |
| Record length in the input | 1000 s |
| Stations | `3Dstations.txt` (81) |

## Mesh

```
python -m salvus_mesh_lite.interface AxiSEM --basic.model Enceladus.bm --basic.period 2 --spherical.min_radius 150 --chunk2D.max_colatitude 60 --output_file enceladus_2s_chunk.e
```

## Run

```
mpirun -np N axisem3d --input input --output output
```

## Files

- `input/`: the five `inparam.*.yaml` files, `3Dstations.txt`, mesh, `Enceladus.bm`, the two model files.
- `make_figures.py`: figures from the station output.
- Fig. 8 as printed (records and maps low-passed at the 2-s mesh period): `../paper_figure_tools/fig8_composite.py`, mode `v6`.
- `paper_figure/make_figures.py`, `pub_style.py`: script of Fig. 8b as submitted. It reads `data/3Dstations.txt` and
  `data/stations/` (copy of the station output; case `ex08` in `../paper_figure_tools/package.sh`, where this
  example is the folder `08_icesheet_regional`). For the paper it was run through
  `../paper_figure_tools/bigfont_runner.py` with `FONT_SCALE=1.5`.

No file is larger than 50 MB. The inputs are those of `examples/08_icesheet_regional/input` of tag `v2.1.0`, unchanged.
