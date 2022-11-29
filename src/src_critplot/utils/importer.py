# -*- coding: utf-8 -*-

import os
from core_gui_atomistic import helpers
from src_critplot.program.topond import atomic_data_from_output
from src_critplot.program.critic2 import atoms_from_xyz, structure_from_cro_file
from src_critplot.utils.critplot_project_file import CritPlotProjectFile


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

            if file_format == "critic_cro":
                models = structure_from_cro_file(filename)
        return models
