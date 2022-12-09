import os
from core_gui_atomistic.gui4dft_project_file import GUI4dftProjectFile
from src_critplot.models.critical_point import CriticalPoint
from src_critplot.models.atomic_model_cp import AtomicModelCP


class CritPlotProjectFile(GUI4dftProjectFile):

    @staticmethod
    def project_file_writer(model):
        text = GUI4dftProjectFile.project_file_writer(model)
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
            atomic_model = GUI4dftProjectFile.project_file_reader(file_name)
            model = AtomicModelCP(atomic_model[0])
            f = open(file_name)
            row = f.readline()
            while row:
                if row.find("%critical points") >= 0:
                    row = f.readline().split()
                    while row.find("%critical points") < 0:
                        if row[0] == "xb":
                            let = "xb"
                            title = let + row[0]
                            x = float(row[1])
                            y = float(row[2])
                            z = float(row[3])
                            charge = -1
                            new_atom = CriticalPoint([x, y, z, let, charge])
                            # new_atom.set_property("title", title)
                            # new_atom.set_property("rho", crit_info[9])
                            # new_atom.set_property("grad", crit_info[10])
                            # new_atom.set_property("lap", crit_info[11])
                            # new_atom.set_property("text", crit_info[12])
                            # new_atom.set_property("atom1", int(data[6]))
                            # new_atom.set_property("atom2", int(data[10]))
                            # translation1 = int(data[7]) * model.lat_vector1 + int(data[8]) * model.lat_vector2 + \
                            #                int(data[9]) * model.lat_vector3
                            # translation2 = int(data[11]) * model.lat_vector1 + int(data[12]) * model.lat_vector2 + \
                            #                int(data[13]) * model.lat_vector3
                            # new_atom.set_property("atom1_translation", translation1)
                            # new_atom.set_property("atom2_translation", translation2)
                            model.add_critical_point(new_atom)
                        print(row)
                        row = f.readline()
                row = f.readline()
            f.close()
        return [model]
