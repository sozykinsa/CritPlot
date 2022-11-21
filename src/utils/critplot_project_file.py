from core_gui_atomistic import helpers
from core_gui_atomistic.gui4dft_project_file import GUI4dftProjectFile
from models.atomic_model_cp import AtomicModelCP


class CritPlotProjectFile(GUI4dftProjectFile):

    @staticmethod
    def project_file_writer(model):
        text = GUI4dftProjectFile.project_file_writer(model)
        try:
            bcp = model.bcp
            text += "%bond critical points\n"
            for point in bcp:
                text += point.to_string() + "\n"
            text += "%end bond critical points\n"
        except:
            pass

        return text

    @staticmethod
    def project_file_reader(file_name):
        atomic_model = GUI4dftProjectFile.project_file_reader(file_name)
        model = AtomicModelCP(atomic_model[0])
        return [model]
