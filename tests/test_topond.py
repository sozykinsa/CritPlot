from program.topond import number_of_atoms_from_outp, get_cell, atomic_data_from_output
import numpy as np
import pytest


def test_number_of_atoms_from_outp(tests_path):
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    n = number_of_atoms_from_outp(f_name)
    assert n == 4


def test_get_cell(tests_path):
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    cell = get_cell(f_name)
    assert len(cell) == 3
    assert np.array(cell[0]) == pytest.approx(np.array([3.57945, 2.34575, 0.00000]))


def test_atomic_data_from_output(tests_path):
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    models = atomic_data_from_output(f_name, True)
    pos = models[0].get_positions()
    assert len(pos) == 9
