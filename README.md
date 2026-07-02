# Example gallery

This gallery holds the **high-resolution meshes** and the **required 3-D model
data** that the examples in `../examples/` reference but that are **too large to
ship in the code repository** (the slim release caps every file at 4 MB and ships
no simulation output). The files here let you run the examples at full resolution.

## How to use it

The gallery mirrors the example layout: data is grouped under one subdirectory
per example, named to match the corresponding `../examples/<name>/` folder.

To run an example at full resolution:

1. Find the example you want under `example_gallery/<example>/`.
2. Copy the file(s) into the **same relative input path** they occupy in the
   matching `../examples/<example>/` folder — the destination is given by the
   **Large data** section of that example's `README.md` (e.g. `input/`,
   `input_share/`, `input_with_3d_crust/`, `input_forward/`, …).
3. Point the example's config at the file if needed: for meshes, set
   `exodus_mesh:` in `inparam.model.yaml` to the mesh you copied in; 3-D model
   files are referenced by their `nc_data_file` / model entries already.
4. Run the solver as described in the example's `README.md`.

Nothing in the gallery is required to *run* an example at its shipped coarse
resolution — each example already ships a working ≤ 4 MB mesh. The gallery is
only needed to reproduce the higher-resolution (headline) runs or to supply the
3-D model data that some configs reference.

## What lives where

| Example | Gallery file(s) | Drop into `examples/<example>/…` | Notes |
|---------|-----------------|----------------------------------|-------|
| `01_S362ANI_EMC_global` | `S362ANI_percent.nc` (~4.9 MB) | `input/` | 3-D mantle model; also on IRIS EMC (model S362ANI). |
| `02_3d_crust_S362ANI_regional` | `S362ANI_percent.nc` (~4.9 MB); `crust1_volumetric.nc` (~30 MB) | `input_share/`; `input_with_3d_crust/` | Crust file also regenerable with `gen_crust1.ipynb`. |
| `04_simple_3d_shapes` | `example_sphere_cartesian.nc` (~70 MB); `plume_1.nc` (~71 MB) | `example_input_cartesian/`; `example_single_plume/` | Both also regenerable in seconds (`gen_blob.py` / `gen_plume.py`). Release-paper model `paper_example.nc` is on Zenodo (4.9 GB), not here. |
| `05_anisotropy_global` | `jons_azi_410_3.nc` (~19 MB); `scanning_output_Nr.nc` (~18 MB); high-res 50 s→finer meshes | `2012-07-03_paper_example_50s/sim_US32_*/input/` | `scanning_output_Nr.nc` only needed by `sim_US32_olivineE` (POINTWISE Nr). |
| `07_kernels` | `S362ANI_radius.nc` (~9.4 MB) | `kernel_3D/input_forward/` | Also regenerable with `convert_s362ani_to_radius.py`. |
| `08_atmosphere_Mars_global` | `tayak_60km_20s.e` (20 s Mars mesh) | `axisem3d_mars_atm/input/` | Also regenerable with `gen_mesh.sh 20`. |
| `09_icesheet_regional` | `volumetric_ice_shell_properties.nc` (~44 MB) | `input/` | Volumetric ice-shell model (Dapré et al., 2024). |
| `10_helioseismology_sun` | `solar_master_epw4.e` (~38 MB, production); `solar_new_test_500.e` (~9.5 MB, quick test) | `input/` | Both derived from `solar_master.bm`; also regenerable with `salvus_mesh_lite`. |
| `11_wave_visualization_Moon` | `moon_issi_m1_4s.e` (~9.8 MB); `moon_issi_m1_2s.e` (~38 MB); `moon_issi_m1_1s.e` (~158 MB) | `sim/input/` | High-resolution lunar meshes; also regenerable by lowering `--basic.period`. |

Examples `00`, `03`, `06`, and `12` do **not** require gallery data (see the
**Zenodo** section for `03`/`06`).

## Very large datasets — Zenodo

A few datasets are too large even for the gallery and are hosted on Zenodo. They
are **not** stored here; download them and place them into the matching example's
input directory as documented in that example's `README.md`:

- **Example 06 (`06_finite_source_local`)** — the SFBA 3-D velocity layers
  `SFBA_layer1_min500.nc`, `SFBA_layer2.nc`, `SFBA_layer3.nc`, `SFBA_layer4.nc`,
  and `topo_smooth_100m_upsampled.nc`. Zenodo: https://zenodo.org/records/20017355
- **Example 03 (`03_salt_body_SEG_local`)** — the SEG/EAGE C3 model archive
  `SEG_C3_data.tar.bz2` (`SEG_C3_SOLID.nc` ~1 GB + `SEG_C3_UNDULATION.nc`); place
  it in `input3D/` and `run3D.sh` extracts it automatically.
- **Example 04 (`04_simple_3d_shapes`)** — the 4.9 GB release-paper scattering
  model `paper_example.nc` (also regenerable with `generate_3D_model.py` from the
  shipped `randomblobs` list).
