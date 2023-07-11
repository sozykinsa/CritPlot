# -*- coding: utf-8 -*-
# Python 3
from typing import List

import math
import numpy as np
from copy import deepcopy
from numpy.linalg import inv
from src_critplot.models.critical_point import CriticalPoint
from core_gui_atomistic.atom import Atom
from core_gui_atomistic.atomic_model import AtomicModel
from core_gui_atomistic import helpers


class AtomicModelCP(AtomicModel):
    def __init__(self, new_atoms: List[Atom] = []):
        super().__init__(new_atoms)
        self.cps: List[CriticalPoint] = []

    def n_cps(self):
        lets = ("xb", "xr", "xc", "nn")
        n = 0
        for cp in self.cps:
            if cp.let in lets:
                n += 1
        return n

    def n_bcp(self):
        # return len(self.cps)
        n = 0
        for cp in self.cps:
            if cp.let == "xb":
                n += 1
        return n

    def n_ccp(self):
        n = 0
        for cp in self.cps:
            if cp.let == "xc":
                n += 1
        return n

    def n_rcp(self):
        n = 0
        for cp in self.cps:
            if cp.let == "xr":
                n += 1
        return n

    def n_ncp(self):
        n = 0
        for cp in self.cps:
            if cp.let == "nn":
                n += 1
        return n

    def bond_path_points_optimize(self) -> None:
        """Remove redundant points from all bond paths."""
        for cp in self.cps:
            if cp.let == "xb":
                bond1 = cp.bonds.get("bond1")
                bond2 = cp.bonds.get("bond2")
                if (bond1 is not None) and (bond2 is not None):
                    self.critical_path_simplifier("bond1", cp)
                    self.critical_path_simplifier("bond2", cp)

    @staticmethod
    def critical_path_simplifier(b, cp) -> None:
        """Remove redundant points on bond path (they lie on the same straight line).
        b: name of bond path (bond1 or bond2)
        cp: critical point
        """
        bond = deepcopy(cp.bonds.get(b))
        if bond is None:
            return
        # i = 2
        # while (i < len(bond)) and (len(bond) > 1):
        #     m = (bond[i].x - bond[i - 1].x) * (bond[i - 2].y - bond[i - 1].y) * (bond[i - 2].z - bond[i - 1].z)
        #     j = (bond[i].y - bond[i - 1].y) * (bond[i - 2].x - bond[i - 1].x) * (bond[i - 2].z - bond[i - 1].z)
        #     k = (bond[i].z - bond[i - 1].z) * (bond[i - 2].x - bond[i - 1].x) * (bond[i - 2].y - bond[i - 1].y)
        #     i += 1
        #     if (math.fabs(m - j) < 1e-6) and (math.fabs(m - k) < 1e-6):
        #         bond.pop(i - 2)
        #         i -= 1

        if len(bond) > 2:
            k = 3
            i = 0
            while (i < len(bond)) and (len(bond) > 1):
                for j in range(1, k):
                    if (i + 1 < len(bond)):
                        bond.pop(i + 1)
                for j in range(k):
                    if (i + 2 < len(bond)):
                        bond.pop(i + 2)
                i += 2

        cp.bonds[b + "opt"] = bond

    def move(self, dl: np.ndarray):
        """Move model by the vector."""
        super().move(dl)
        for cp in self.cps:
            cp.xyz += dl
            self.move_bond_path(dl, cp.bonds.get("bond1"))
            self.move_bond_path(dl, cp.bonds.get("bond2"))
            self.move_bond_path(dl, cp.bonds.get("bond1opt"))
            self.move_bond_path(dl, cp.bonds.get("bond2opt"))
        return self.atoms

    @staticmethod
    def move_bond_path(dl, bond):
        if bond:
            for bp in bond:
                bp.xyz += dl

    def go_to_positive_coordinates_translate(self):
        self.translated_atoms_remove()
        self.go_to_positive_array_translate(self.atoms)
        self.find_bonds_fast()
        self.go_to_positive_array_translate(self.cps)
        self.bond_path_opt_update()

    def bond_path_opt_update(self):
        for i in range(len(self.cps)):
            if self.cps[i].let == "xb":
                if "bond1" not in self.cps[i].bonds:
                    print("strange bcp: ", self.cps[i])
                else:
                    self.cps[i].bonds.pop("bond1")
                    self.cps[i].bonds.pop("bond2")
                    self.cps[i].bonds.pop("bond1opt")
                    self.cps[i].bonds.pop("bond2opt")

                    ind1 = self.cps[i].get_property("atom1")
                    ind2 = self.cps[i].get_property("atom2")
                    if (ind1 <= len(self.atoms)) and (ind2 <= len(self.atoms)):
                        p1 = CriticalPoint([*self.atoms[ind1 - 1].xyz, "xz", 1])
                        p2 = CriticalPoint([*self.cps[i].xyz, "xz", 1])
                        p3 = CriticalPoint([*self.atoms[ind2 - 1].xyz, "xz", 1])

                        self.add_bond_path_point([p2, p1])
                        self.add_bond_path_point([p2, p3])
        self.bond_path_points_optimize()

    def move_atoms_to_cell(self):
        a_inv = inv(self.lat_vectors)
        self.move_object_to_cell(self.atoms, a_inv)
        self.move_object_to_cell(self.cps, a_inv)
        for cp in self.cps:
            bond1 = cp.bonds.get("bond1")
            bond2 = cp.bonds.get("bond2")
            if bond1 is not None and bond2 is not None:
                bond1opt = cp.bonds.get("bond1opt")
                bond2opt = cp.bonds.get("bond2opt")
                self.move_object_to_cell(bond1, a_inv)
                self.move_object_to_cell(bond2, a_inv)
                self.move_object_to_cell(bond1opt, a_inv)
                self.move_object_to_cell(bond2opt, a_inv)

    def go_to_positive_coordinates(self):
        dr = np.array([self.minX(), self.minY(), self.minZ()])
        self.go_to_positive_array(self.atoms, dr)
        self.go_to_positive_array(self.cps, dr)
        for cp in self.cps:
            bond1 = cp.bonds.get("bond1")
            bond2 = cp.bonds.get("bond2")
            if bond1 is not None and bond2 is not None:
                bond1opt = cp.bonds.get("bond1opt")
                bond2opt = cp.bonds.get("bond2opt")
                self.go_to_positive_array(bond1, dr)
                self.go_to_positive_array(bond2, dr)
                self.go_to_positive_array(bond1opt, dr)
                self.go_to_positive_array(bond2opt, dr)

    def convert_from_direct_to_cart(self):
        super().convert_from_direct_to_cart()
        for cp in self.cps:
            cp.xyz = np.dot(self.lat_vectors, cp.xyz)

    def add_critical_point(self, cp):
        self.cps.append(deepcopy(cp))

    def delete_all_bond_paths(self):
        for cp in self.cps:
            cp.bonds["bond1"] = None
            cp.bonds["bond2"] = None

    def add_bond_path_point(self, points):
        for cp in self.cps:
            distance = math.dist(cp.xyz, points[0].xyz)
            if distance < 1e-4:
                if cp.bonds.get("bond1") is None:
                    cp.bonds["bond1"] = deepcopy(points)
                else:
                    cp.bonds["bond2"] = deepcopy(points)

    @staticmethod
    def atoms_of_bond_path(cp):
        ind1 = cp.get_property("atom1")
        ind2 = cp.get_property("atom2")
        return ind1, ind2

    def create_csv_file_cp(self, cp_list, show_bcp, show_ccp, show_rcp, show_ncp, show_nnatr, delimiter: str = ";"):
        text = ""

        if show_bcp:
            cp_type = "xb"
            data, title = self.csv_for_cp_type(cp_list, cp_type, delimiter)
            text += title + "\n" + data + "\n"

        if show_ccp:
            cp_type = "xc"
            data, title = self.csv_for_cp_type(cp_list, cp_type, delimiter)
            text += title + "\n" + data + "\n"

        if show_rcp:
            cp_type = "xr"
            data, title = self.csv_for_cp_type(cp_list, cp_type, delimiter)
            text += title + "\n" + data + "\n"

        if show_nnatr:
            cp_type = "nn"
            data, title = self.csv_for_cp_type(cp_list, cp_type, delimiter)
            text += title + "\n" + data + "\n"

        if show_ncp:
            cp_type = "nucleus"
            data, title = self.csv_for_cp_type(cp_list, cp_type, delimiter)
            text += title + "\n" + data + "\n"
        return text

    def csv_for_cp_type(self, cp_list, cp_type, delimiter):
        data = ""
        title = ""
        # print(cp_type)
        for ind in cp_list:
            cp = self.cps[ind - 1]
            if (cp.let == cp_type) or (cp_type == "nucleus") and (cp.let not in ["xb", "xc", "xr", "nn"]):
                title = "CP" + delimiter
                data += str(ind) + delimiter
                prop_list = ["rho", "grad", "lap"]
                for prop in prop_list:
                    title += prop + delimiter
                    data += cp.get_property(prop) + delimiter

                title += "atoms" + delimiter + "dist" + delimiter
                atom_to_atom = cp.get_property("atom_to_atom")
                if atom_to_atom is None:
                    data += '" ' + "-" + ' "' + delimiter + '" ' + "-" + ' "' + delimiter
                else:
                    data += atom_to_atom + delimiter
                    cp_bp_len = cp.get_property("cp_bp_len")
                    dist_line = round(cp_bp_len, 4)
                    data += '" ' + str(dist_line) + ' "' + delimiter
                data_list = cp.get_property("text")
                if data_list:
                    data = self.text_field_to_csv(data, data_list, delimiter, title)
                    title = self.text_field_to_csv_title(data_list, delimiter, title)
        return data, title

    @staticmethod
    def text_field_to_csv(data, data_list, delimiter, title):
        data_list = data_list.split("\n")
        i = 0
        while i < len(data_list):
            if (data_list[i].find("Hessian") < 0) and (len(data_list[i]) > 0):
                col_data = data_list[i].split(":")[1].split()
                for k in range(len(col_data)):
                    data += '"' + col_data[k] + '"' + delimiter
                i += 1
            else:
                i += 4
        data += '\n'
        return data

    @staticmethod
    def text_field_to_csv_title(data_list, delimiter, title):
        data_list = data_list.split("\n")
        i = 0
        while i < len(data_list):
            if (data_list[i].find("Hessian") < 0) and (len(data_list[i]) > 0):
                col_title = helpers.spacedel(data_list[i].split(":")[0])
                col_data = data_list[i].split(":")[1].split()
                for k in range(len(col_data)):
                    title += '"' + col_title
                    if len(col_data) > 1:
                        title += "_" + str(k + 1)
                    title += '"' + delimiter
                i += 1
            else:
                i += 4
        return title
