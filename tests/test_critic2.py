from utils.import_export import ImporterExporter
from models.atomic_model_cp import AtomicModelCP


def test_atoms_of_bond_path(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = ImporterExporter.import_from_file(f_name)
    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[1])
    assert atom1 is None
    assert atom2 is None

    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[2])
    assert atom1 is None
    assert atom2 is None

    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[3])
    assert atom1 == 3
    assert atom2 == 3

    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[4])
    assert atom1 == 3
    assert atom2 == 3

    atom1, atom2 = AtomicModelCP.atoms_of_bond_path(model[0].cps[5])
    assert atom1 == 2
    assert atom2 == 1
