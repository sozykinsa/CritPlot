from core_atomistic.project_file import ProjectFile


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
