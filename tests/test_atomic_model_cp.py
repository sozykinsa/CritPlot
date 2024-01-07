from models.cp_model import AtomicModelCP
from programs.topond import TopondModelCP
from programs.critic2 import Critic2ModelCP, parse_bondpaths


def test_critical_path_simplifier(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model: AtomicModelCP = Critic2ModelCP(f_name)
    assert len(model.cps) == 11
    bond1 = model.cps[7].bonds.get("bond1")
    assert len(bond1) == 2
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    models = parse_bondpaths(f_name, model)
    assert len(models[0].cps) == 11
    assert models[0].cps[7].get_property("atom1") == 2
    bond1 = models[0].cps[7].bonds.get("bond1")
    assert len(bond1) == 44


def test_get_cp_types(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model: AtomicModelCP = Critic2ModelCP(f_name)
    cp_types = model.get_cp_types()
    assert len(cp_types) == 4


def test_create_csv_file_cp(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model: AtomicModelCP = Critic2ModelCP(f_name)
    assert len(model.cps) == 11
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    models = parse_bondpaths(f_name, model)
    cp_list = range(len(models[0].cps))
    text = models[0].create_csv_file_cp(cp_list, True, False, False, False, False)
    assert len(text) == 2086


def test_go_to_positive_coordinates_translate(tests_path):
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    model = TopondModelCP(f_name, True)
    assert len(model.atoms) == 9
    model.go_to_positive_coordinates_translate()
    assert len(model.atoms) == 4
