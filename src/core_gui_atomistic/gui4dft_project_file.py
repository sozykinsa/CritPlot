import os
from core_gui_atomistic import helpers
from core_gui_atomistic.atomic_model import AtomicModel


class GUI4dftProjectFile(object):

    @staticmethod
    def project_file_writer(model):
        text = ""
        try:
            vectors = model.lat_vectors
            text += "%vectors\n"
            for vect in vectors:
                text += str(vect[0]) + " " + str(vect[1]) + " " + str(vect[2]) + "\n"
            text += "%end vectors\n"
        except:
            pass
        try:
            atoms = model.atoms
            text += "%atoms\n"
            text += str(len(atoms)) + "\n\n"
            for point in atoms:
                text += point.to_string() + "\n"
            text += "%end atoms\n"
        except:
            pass

        return text

    @staticmethod
    def project_file_reader(file_name):
        model = AtomicModel()
        if os.path.exists(file_name):
            f = open(file_name)
            row = f.readline()
            while row:
                if row.find("%atoms") >= 0:
                    number_of_atoms = int(f.readline())
                    #row = f.readline()
                    model = AtomicModel.atoms_from_xyz_structure(number_of_atoms, f)
                row = f.readline()
            f.close()
        return [model]

