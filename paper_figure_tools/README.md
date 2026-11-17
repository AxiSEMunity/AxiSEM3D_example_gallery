# Shared figure tools

- `fig8_composite.py`: Fig. 8 as printed (`python fig8_composite.py v6 out.pdf`): (a) Mars, (b) Enceladus, (c) the Sun,
  with every record and snapshot low-passed at the mesh period of its run (`fig8_lowpass.py`). It reads the three runs
  from one folder (`FIG8_ROOT`): `mars/`, `enceladus/`, `sun/` (runs of the input folders of Examples 8, 9 and 10, each
  with its output in `output/`) and
  `sun_xs.npz`, written by `dump_sun_v6.py` from `sun/output/elements/sun_cross_section`. Details in its docstring.
- `pdf_plain_images.py`: rewrites the images of a figure PDF as plain opaque RGB images (Matplotlib writes palette
  images, which some PDF conversions show as black boxes). Used for Figs. 4 to 8.
- `package.sh`, `bigfont_runner.py`: extraction of compact data files and the large-font run of the figure scripts
  as submitted.
