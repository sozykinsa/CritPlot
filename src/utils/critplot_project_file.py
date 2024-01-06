import os
from src_core_atomistic.project_file import ProjectFile
from models.critical_point import CriticalPoint
from models.atomic_model_cp import AtomicModelCP


class CritPlotProjectFile(ProjectFile):

    @staticmethod
    def project_file_writer(model):
        text = ProjectFile.project_file_writer(model)
        try:
            cps = model.cps
            text += "%critical points\n"
            for point in cps:
                text += point.to_string() + "\n"
            text += "%end critical points\n"
        except:
            pass

        return text

    @staticmethod
    def project_file_reader(file_name):
        model = None
        if os.path.exists(file_name):
            atomic_model = ProjectFile.project_file_reader(file_name)
            model = AtomicModelCP(atomic_model[0])
            f = open(file_name)
            row = f.readline()
            while row:
                if row.find("%critical points") >= 0:
                    while row.find("%critical points") < 0:
                        row = f.readline().split()
                        if row[0] == "xb":
                            let = "xb"
                            title = let + row[0]
                            x = float(row[1])
                            y = float(row[2])
                            z = float(row[3])
                            charge = -1
                            new_atom = CriticalPoint([x, y, z, let, charge])
                            model.add_critical_point(new_atom)
                        row = f.readline()
                row = f.readline()
            f.close()
        return [model]
