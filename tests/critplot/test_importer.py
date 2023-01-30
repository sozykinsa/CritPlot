from src.src_critplot.utils.import_export import ImporterExporter
from core_gui_atomistic import helpers


def test_check_format(tests_path):
    assert helpers.check_format(str(tests_path / 'ref_data' / 'swcnt(8,0)' / "siesta1.out")) == "unknown"
    assert helpers.check_format(str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")) == "critic_cro"


def test_importer_xyz(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = ImporterExporter.import_from_file(f_name)
    print(model)
    for at in model[0].atoms:
        print(at.to_string())
    assert len(model[0].atoms) == 6
    model[0].translated_atoms_remove()
    assert len(model[0].atoms) == 3
