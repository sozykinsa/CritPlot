# -*- coding: utf-8 -*-

import os
from copy import deepcopy
from core_atomistic import helpers
from programs.topond import atomic_data_from_output
from programs.critic2 import structure_from_cro_file
from programs.cp_project_file import CritPlotProjectFile
from models.cp_model import AtomicModelCP


class ImporterExporter(object):

    @staticmethod
    def import_from_file(filename: str, is_add_translations = False):
        """import file"""
        models = []
        is_critic_open = False
        if os.path.exists(filename):
            file_format = helpers.check_format(filename)
            print("format: ", file_format)
            if file_format == "topond_out":
                models = atomic_data_from_output(filename, is_add_translations)
            elif file_format == "critic_cro":
                models = structure_from_cro_file(filename)
                is_critic_open = True
            elif file_format == "project":
                models = CritPlotProjectFile.project_file_reader(filename)
            else:
                print("Wrong format")
        return models, is_critic_open

    @staticmethod
    def model_for_export(model: AtomicModelCP):
        new_model = deepcopy(model)
        new_model.translated_atoms_remove()
        return new_model

    @staticmethod
    def export_to_file(model, file_name):  # pragma: no cover
        new_model = ImporterExporter.model_for_export(model)
        if file_name.endswith(".data"):
            text = CritPlotProjectFile.project_file_writer(new_model)
            helpers.write_text_to_file(file_name, text)
