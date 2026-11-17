"""Run one of the per-body overview scripts with all font sizes scaled by FONT_SCALE and,
at save time, the same clean-up that was done by hand in Illustrator for the paper's
combined figure: no suptitle, no per-panel title, record-section title only on Mars
('Record section'), colourbar only on Mars."""
import os, sys, runpy
SCALE = float(os.environ.get('FONT_SCALE', '1.5'))
import matplotlib; matplotlib.use('Agg')
import matplotlib.text as mtext, matplotlib.figure as mfig
_orig = mtext.Text.set_fontsize
def _sf(self, fontsize):
    try: fontsize = float(fontsize) * SCALE
    except (TypeError, ValueError): pass
    return _orig(self, fontsize)
mtext.Text.set_fontsize = _sf; mtext.Text.set_size = _sf
_osave = mfig.Figure.savefig
def _save(self, fname, *a, **k):
    name = os.path.basename(str(fname))
    if any(t in name for t in ('mars_atm_section', 'ex08_enceladus_overview', 'sun_overview')):
        mars = 'mars' in name
        if getattr(self, '_suptitle', None) is not None: self._suptitle.set_visible(False)
        for ax in self.axes:
            ttl = ax.get_title()
            if 'record section' in ttl.lower(): ax.set_title('Record section' if mars else '')
            pos = ax.get_position()
            is_cbar = ax.get_label() == '<colorbar>' or (pos.height < 0.045 and pos.width > 0.15) or (pos.width < 0.03 and pos.height > 0.2)
            if is_cbar and not mars: ax.set_visible(False)
    return _osave(self, fname, *a, **k)
mfig.Figure.savefig = _save
script = sys.argv[1]; _d = os.path.dirname(os.path.abspath(script)); os.chdir(_d); sys.path.insert(0, _d); sys.argv = [script]
runpy.run_path(script, run_name='__main__')
