# -*- coding: utf-8 -*-

import os
from core_gui_atomistic import helpers
from src_critplot.program.topond import atomic_data_from_output
from src_critplot.program.critic2 import structure_from_cro_file
from src_critplot.utils.critplot_project_file import CritPlotProjectFile


class Importer(object):

    @staticmethod
    def import_from_file(filename):
        """import file"""
        models = []
        is_critic_open = False
        if os.path.exists(filename):
            file_format = helpers.check_format(filename)
            if file_format == "topond_out":
                models = atomic_data_from_output(filename)
            elif file_format == "critic_cro":
                models = structure_from_cro_file(filename)
                is_critic_open = True
            elif file_format == "project":
                models = CritPlotProjectFile.project_file_reader(filename)
            else:
                print("Wrong format")
        return models, is_critic_open
