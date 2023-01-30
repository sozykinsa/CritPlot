from pathlib import Path
from core_gui_atomistic.atomic_model import AtomicModel
from core_gui_atomistic.periodic_table import TPeriodTable
from src.src_critplot.qtbased.pyqtgraphwidget import PyqtGraphWidget
from src.src_critplot.qtbased.pyqtgraphwidgetimage import PyqtGraphWidgetImage
from src.src_critplot.models.atomic_model_cp import AtomicModelCP
from qtpy.QtWidgets import QApplication


import pytest


@pytest.fixture
def qsapp(qapp_session) -> QApplication:
    yield qapp_session
    qapp_session.processEvents()
    qapp_session.closeAllWindows()


@pytest.fixture(scope='session')
def qapp_session(qapp):
    return qapp


@pytest.fixture
def tests_path() -> Path:
    return Path(__file__).parent


@pytest.fixture
def h2o_model() -> AtomicModel:
    #         x    y    z   let charge
    atoms = [[0.0, 0.0, 0.0, "O", 8], [-1.0, 0.0, 0.0, "H", 1], [1.0, 0.0, 0.0, "H", 1]]
    return AtomicModel(atoms)


@pytest.fixture
def h2o_model_cp() -> AtomicModelCP:
    #         x    y    z   let charge
    atoms = [[0.0, 0.0, 0.0, "O", 8], [-1.0, 0.0, 0.0, "H", 1], [1.0, 0.0, 0.0, "H", 1]]
    return AtomicModelCP(atoms)


@pytest.fixture
def swnt_33() -> Path:
    #         x    y    z   let charge
    atoms = [[-0.70056811, 1.92479506, - 0.61920816, "C", 6], [-1.31663736, - 1.56910731, - 0.61920816, "C", 6],
             [-2.01720547, - 0.35568775, - 0.61920816, "C", 6], [2.01720547, - 0.35568775, - 0.61920816, "C", 6],
             [1.31663736, - 1.56910731, - 0.61920816, "C", 6], [0.70056811, 1.92479506, - 0.61920816, "C", 6],
             [-1.31663736, 1.56910731, 0.61920816, "C", 6], [-0.70056811, - 1.92479506, 0.61920816, "C", 6],
             [-2.01720547, 0.35568775, 0.61920816, "C", 6], [2.01720547, 0.35568775, 0.61920816, "C", 6],
             [0.70056811, - 1.92479506, 0.61920816, "C", 6], [1.31663736, 1.56910731, 0.61920816, "C", 6]]
    return AtomicModel(atoms)


@pytest.fixture
def period_table() -> TPeriodTable:
    return TPeriodTable()


@pytest.fixture
def get_graph_widget(qsapp):

    def factory_function():
        widget = PyqtGraphWidget()
        widget.show()
        return widget

    return factory_function


@pytest.fixture
def graph_widget(get_graph_widget):
    return get_graph_widget()


@pytest.fixture
def get_graph_image_widget(qsapp):

    def factory_function():
        widget = PyqtGraphWidgetImage()
        widget.show()
        return widget

    return factory_function


@pytest.fixture
def graph_image_widget(get_graph_image_widget):
    return get_graph_image_widget()
