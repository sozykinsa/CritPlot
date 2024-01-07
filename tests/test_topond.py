from programs.topond import TopondModelCP
import numpy as np
import pytest


def test_topond_model_cp_constructor(tests_path):
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    model = TopondModelCP(f_name, True)
    n = model.n_atoms()
    assert n == 4

    cell = model.get_cell(f_name)
    assert len(cell) == 3
    assert np.array(cell[0]) == pytest.approx(np.array([3.57945, 2.34575, 0.00000]))

    pos = model.get_positions()
    assert len(pos) == 9
