# -*- coding: utf-8 -*-

import os
from program.topond import atomic_data_from_output
from program.critic2 import atoms_from_xyz
from core_gui_atomistic import helpers
from utils.critplot_project_file import CritPlotProjectFile


class Importer(object):

    @staticmethod
    def import_from_file(filename):
        """import file"""
        models = []
        if os.path.exists(filename):
            file_format = helpers.check_format(filename)

            if file_format == "SiestaXYZ":
                models = atoms_from_xyz(filename)

            if file_format == "topond_out":
                models = atomic_data_from_output(filename)

            if file_format == "project":
                models = CritPlotProjectFile.project_file_reader(filename)
        return models
