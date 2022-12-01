import pytest
from src.src_critplot.program.critic2 import parse_cp_properties


def test_critic2_section(critplot_application, tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    window = critplot_application
    window.menu_open(f_name)
    model = window.models[-1]
    assert len(model.bcp) == 5
    assert "rho" not in model.bcp[0].properties
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    parse_cp_properties(f_name, model)
    assert float(model.bcp[0].properties["rho"]) == pytest.approx(1.68288239)
    window.add_cp_to_list()
    assert window.ui.FormCPlist.count() == 0
    window.selected_cp_changed(1)
    assert window.ui.selectedCP_nuclei.text() == "O3-O3"
    window.ui.selectedCP.setText("2")
    window.add_cp_to_list()
    window.add_cp_to_list()
    assert window.ui.FormCPlist.item(0).text() == "2"
    window.ui.FormCPlist.setCurrentRow(0)
    assert window.ui.FormCPlist.count() == 1
    window.delete_cp_from_list()
    assert window.ui.FormCPlist.count() == 0


def test_create_cri_file(critplot_application, tests_path):
    f_name = str(tests_path / 'ref_data' / 'critic2' / "cp-file.xyz")
    window = critplot_application
    window.menu_open(f_name)
    model = window.models[-1]
    f_name = str(tests_path / 'ref_data' / 'critic2' / "siesta-1-cp.cro")
    parse_cp_properties(f_name, model)

    window.ui.selectedCP.setText("1")
    window.add_cp_to_list()
    assert len(window.models[-1].bcp) == 5
    window.delete_cp_from_model()
    assert len(window.models[-1].bcp) == 4
    window.ui.selectedCP.setText("2")
    window.add_cp_to_list()
    window.leave_cp_in_model()
    assert len(window.models[-1].bcp) == 1
