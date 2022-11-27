from src.src_critplot.utils.importer import Importer


def test_importer_xyz(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    model = Importer.import_from_file(f_name)
    assert len(model[0].atoms) == 3
