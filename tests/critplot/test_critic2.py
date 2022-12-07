from src.src_critplot.utils.importer import Importer
from src.src_critplot.models.atomic_model_cp import AtomicModelCP


"""
def test_create_critic2_xyz_file(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model = Importer.import_from_file(f_name)
    assert len(model[0].cps) == 5
    text = create_critic2_xyz_file(model[0].cps, (model[0].cps[1], model[0].cps[2]), True, model[0])
    assert len(text) == 3739
    text = create_critic2_xyz_file(model[0].cps, (model[0].cps[1], model[0].cps[2]), False, model[0])
    assert len(text) == 8375
"""


def test_atoms_of_bond_path(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = Importer.import_from_file(f_name)
    # for cp in model[0].cps:
    #     print(cp.let, cp.xyz, cp.bonds, cp.properties)
    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[1])
    assert atom1 == 3
    assert atom2 == 3

    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[2])
    assert atom1 == 2
    assert atom2 == 1


"""
def test_create_csv_file_cp(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = Importer.import_from_file(f_name)
    text = model[0].create_csv_file_cp([0, 1], delimiter=";")
    assert len(text) > 4


def test_check_cro_file(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    box_bohr, box_ang, box_deg, cps = check_cro_file(f_name)
    assert len(box_bohr) == 3
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp-error.cro")
    box_bohr, box_ang, box_deg, cps = check_cro_file(f_name)
    assert len(box_bohr) == 0


def test_create_cri_file(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = Importer.import_from_file(f_name)
    textl, lines, te, text = create_cri_file((1, 2), 1, True, model[0], "")
    assert len(textl) == 234
    textl, lines, te, text = create_cri_file((1, 2), 1, False, model[0], "")
    assert len(textl) == 234
"""
