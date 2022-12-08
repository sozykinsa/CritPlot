from src.src_critplot.utils.import_export import ImporterExporter


def test_importer_xyz(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    model, fl = ImporterExporter.import_from_file(f_name)
    assert len(model[0].atoms) == 3
