from utils.importer import Importer


def test_importer_xyz(tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    model = Importer.import_from_file(f_name)
    assert len(model[0].atoms) == 3


#def test_importer_outp(tests_path):
#    f_name = str(tests_path / 'ref_data' / 'topond' / "topond-I.outp")
#    model = Importer.import_from_file(f_name)
#    assert len(model[0].atoms) == 1
