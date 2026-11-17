"""axikernels ElementOutput for element output with ALL GLL points (GLL_points_one_edge: FULL).

The published reader supports only GLL_points_one_edge = [0, 2, 4] (a 3 x 3 subset of the 5 x 5 nodes of each element)
and interpolates with quadratic Lagrange polynomials, while the solution is a degree-4 polynomial in each element: this
leaves element-scale ripples in kernels. Here the element geometry (point location) still uses the 3 x 3 subset of the
stored coordinates, unchanged, and the DATA are interpolated with the full degree-4 Lagrange basis on the solver's nodes:
GLL (-1, -sqrt(3/7), 0, sqrt(3/7), 1) and, on axial elements in the s direction, GLJ (AxiSEM3D spectrals.hpp, NPOL 4).
Node order ipnt = ipol * 5 + jpol (ipol along xi, jpol along eta), as for the 3 x 3 output.
Use: from axk_full import ElementOutputFull as ElementOutput (same interface)."""
import numpy as np
from axikernels.core.handlers import element_output as _eo

GLL5 = np.array([-1.0, -np.sqrt(3.0 / 7.0), 0.0, np.sqrt(3.0 / 7.0), 1.0])
GLJ5 = np.array([-1.0, -0.507787629558, 0.132300820777, 0.708820142114, 1.0])
SUB9 = np.array([5 * i + j for i in (0, 2, 4) for j in (0, 2, 4)])
_orig_weights = _eo.interpolation_weights_9node


def _lagrange(t, nodes):
    w = np.ones(len(nodes))
    for i, xi in enumerate(nodes):
        for j, xj in enumerate(nodes):
            if j != i:
                w[i] *= (t - xj) / (xi - xj)
    return w


def _weights25(xi, eta, xi_nodes, eta_nodes):
    """Called by the reader with the 3-node subsets; the middle xi node tells GLL (0) from GLJ (0.1323)."""
    xn = GLJ5 if abs(float(xi_nodes[1])) > 1e-6 else GLL5
    return np.outer(_lagrange(xi, xn), _lagrange(eta, GLL5)).ravel()


class ElementOutputFull(_eo.ElementOutput):
    def __init__(self, path_to_element_output):
        super().__init__(path_to_element_output=path_to_element_output)
        self._full = {}
        for g in self.element_groups:
            md = self.element_groups_info[g]['metadata']; inp = self.element_groups_info[g]['inplane']
            c = np.asarray(md['list_element_coords'])
            full = c.ndim == 3 and c.shape[1] == 25
            self._full[g] = full
            if full:
                md['list_element_coords_full'] = c
                md['list_element_coords'] = c[:, SUB9, :]
                inp['GLL_points_one_edge_stored'] = inp['GLL_points_one_edge']
                inp['GLL_points_one_edge'] = [0, 2, 4]          # passes the reader's check; data use all 25 nodes
        print('ElementOutputFull:', {g: ('25 nodes' if f else '9 nodes') for g, f in self._full.items()}, flush=True)

    def load_data_from_element_group(self, points, group, channel_slices=None, time_slices=None, pbar=None):
        if not self._full.get(group, False):
            return super().load_data_from_element_group(points, group, channel_slices, time_slices, pbar)
        _eo.interpolation_weights_9node = _weights25
        try:
            return super().load_data_from_element_group(points, group, channel_slices, time_slices, pbar)
        finally:
            _eo.interpolation_weights_9node = _orig_weights
