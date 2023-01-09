from src_critplot.models.atomic_model_cp import AtomicModelCP
from src_critplot.program.topond import atomic_data_from_output
from src_critplot.program.critic2 import structure_from_cro_file, parse_bondpaths


def test_critical_points(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    models: AtomicModelCP = structure_from_cro_file(f_name)
    assert len(models[0].cps) == 11
    assert models[0].bond_path_len(models[0].cps[5]) == 3.5165430594833227
    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
    models = atomic_data_from_output(f_name)
    assert len(models[0].cps) == 7
    assert models[0].bond_path_len(models[0].cps[3]) == 4.4238456974571685


def test_critical_path_simplifier(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    models: AtomicModelCP = structure_from_cro_file(f_name)
    assert len(models[0].cps) == 11
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    models = parse_bondpaths(f_name, models[0])
    assert len(models[0].cps) == 11
    assert models[0].cps[7].get_property("atom1") == 2
    bond1 = models[0].cps[7].bonds.get("bond1")
    assert len(bond1) == 44
