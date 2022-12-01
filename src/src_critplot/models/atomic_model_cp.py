# -*- coding: utf-8 -*-
# Python 3

import math
import numpy as np
from copy import deepcopy
from numpy.linalg import inv
from core_gui_atomistic.atom import Atom
from core_gui_atomistic.atomic_model import AtomicModel
from core_gui_atomistic import helpers


class AtomicModelCP(AtomicModel):
    def __init__(self, new_atoms: list = []):
        super().__init__(new_atoms)
        self.cps = []

    def bond_path_points_optimize(self):
        i = 0

        while i < len(self.cps):
            bond1 = self.cps[i].get_property("bond1")
            bond2 = self.cps[i].get_property("bond2")
            if (bond1 is None) or (bond2 is None):
                self.cps.pop(i)
                i -= 1
            else:
                if (bond1[-1].x == bond2[-1].x) and (bond1[-1].y == bond2[-1].y) and (bond1[-1].z == bond2[-1].z):
                    self.cps.pop(i)
                    i -= 1
            i += 1

        for cp in self.cps:
            self.critical_path_simplifier("bond1", cp)
            self.critical_path_simplifier("bond2", cp)

    @staticmethod
    def critical_path_simplifier(b, cp):
        bond = deepcopy(cp.get_property(b))
        if bond is None:
            return
        i = 2
        while (i < len(bond)) and (len(bond) > 1):
            m = (bond[i].x - bond[i - 1].x) * (bond[i - 2].y - bond[i - 1].y) * (bond[i - 2].z - bond[i - 1].z)
            j = (bond[i].y - bond[i - 1].y) * (bond[i - 2].x - bond[i - 1].x) * (bond[i - 2].z - bond[i - 1].z)
            k = (bond[i].z - bond[i - 1].z) * (bond[i - 2].x - bond[i - 1].x) * (bond[i - 2].y - bond[i - 1].y)
            i += 1
            if (math.fabs(m - j) < 1e-6) and (math.fabs(m - k) < 1e-6):
                bond.pop(i - 2)
                i -= 1
        cp.set_property(b + "opt", bond)

    def move(self, l_x, l_y, l_z):
        """Move model by the vector."""
        super().move(l_x, l_y, l_z)
        dl = np.array([l_x, l_y, l_z])

        for cp in self.cps:
            cp.xyz += dl
            self.move_bond_path(dl, cp.get_property("bond1"))
            self.move_bond_path(dl, cp.get_property("bond2"))
            self.move_bond_path(dl, cp.get_property("bond1opt"))
            self.move_bond_path(dl, cp.get_property("bond2opt"))
        return self.atoms

    @staticmethod
    def move_bond_path(dl, bond):
        if bond:
            for bp in bond:
                bp.xyz += dl

    def go_to_positive_coordinates_translate(self):
        self.go_to_positive_array_translate(self.atoms)
        self.go_to_positive_array_translate(self.cps)
        self.go_to_positive_array_translate(self.rcp)
        self.go_to_positive_array_translate(self.ccp)
        self.bond_path_opt_update()

    def bond_path_opt_update(self):
        for i in range(len(self.cps)):
            self.cps[i].properties.pop("bond1")
            self.cps[i].properties.pop("bond2")
            self.cps[i].properties.pop("bond1opt")
            self.cps[i].properties.pop("bond2opt")

            ind1 = self.cps[i].get_property("atom1")
            ind2 = self.cps[i].get_property("atom2")
            p1 = Atom([*self.atoms[ind1].xyz, "xz", 1])
            p2 = Atom([*self.cps[i].xyz, "xz", 1])
            p3 = Atom([*self.atoms[ind2].xyz, "xz", 1])

            self.add_bond_path_point([p2, p1])
            self.add_bond_path_point([p2, p3])
        self.bond_path_points_optimize()

    def n_bcp(self):
        return len(self.cps)

    def move_atoms_to_cell(self):
        a_inv = inv(self.lat_vectors)
        self.move_object_to_cell(self.atoms, a_inv)
        self.move_object_to_cell(self.cps, a_inv)
        self.move_object_to_cell(self.rcp, a_inv)
        self.move_object_to_cell(self.ccp, a_inv)

    def go_to_positive_coordinates(self):
        xm = self.minX()
        ym = self.minY()
        zm = self.minZ()
        self.go_to_positive_array(self.atoms, xm, ym, zm)
        self.go_to_positive_array(self.cps, xm, ym, zm)

    def convert_from_direct_to_cart(self):
        super().convert_from_direct_to_cart()
        for cp in self.cps:
            cp.xyz = np.dot(self.lat_vectors, cp.xyz)

    def add_critical_point(self, atom):
        self.cps.append(deepcopy(atom))

    def add_bond_path_point(self, points):
        for cp in self.cps:
            distance = math.dist(cp.xyz, points[0].xyz)
            if distance < 1e-4:
                if cp.get_property("bond1") is None:
                    cp.set_property("bond1", deepcopy(points))
                else:
                    cp.set_property("bond2", deepcopy(points))
                    ind1, ind2 = self.atoms_of_bond_path(self.atoms, self.cps[self.cps.index(cp)])
                    cp.set_property("atom1", ind1)
                    cp.set_property("atom2", ind2)

    @staticmethod
    def atoms_of_bond_path(atoms, cp):
        bond1 = cp.get_property("bond1")
        bond2 = cp.get_property("bond2")
        if bond1 is None or bond2 is None:
            return 0, 0
        minr1 = math.dist(bond1[-1].xyz, atoms[0].xyz)
        minr2 = math.dist(bond2[-1].xyz, atoms[0].xyz)
        ind1 = 0
        ind2 = 0
        for i in range(0, len(atoms)):
            d1 = math.dist(bond1[-1].xyz, atoms[i].xyz)
            d2 = math.dist(bond2[-1].xyz, atoms[i].xyz)

            if d1 < minr1:
                minr1 = d1
                ind1 = i

            if d2 < minr2:
                minr2 = d2
                ind2 = i
        return ind1, ind2

    def create_csv_file_cp(self, cp_list, delimiter: str = ";"):
        title = ""
        data = ""
        for ind in cp_list:
            title = ""
            cp = self.cps[ind]
            title += "BCP" + delimiter + "atoms" + delimiter + "dist" + delimiter
            data += str(ind) + delimiter
            ind1 = cp.get_property("atom1")
            ind2 = cp.get_property("atom2")
            atom1 = self.atoms[ind1].let + str(ind1 + 1)
            atom2 = self.atoms[ind2].let + str(ind2 + 1)
            data += atom1 + "-" + atom2 + delimiter

            dist_line = round(self.atom_atom_distance(ind1, ind2), 4)
            data += '" ' + str(dist_line) + ' "' + delimiter

            data_list = cp.get_property("text")
            if data_list:
                data_list = data_list.split("\n")
                i = 0

                while i < len(data_list):
                    if (data_list[i].find("Hessian") < 0) and (len(data_list[i]) > 0):
                        col_title = helpers.spacedel(data_list[i].split(":")[0])
                        col_data = data_list[i].split(":")[1].split()
                        for k in range(len(col_data)):
                            title += col_title
                            if len(col_data) > 1:
                                title += "_" + str(k + 1)
                            title += delimiter
                            data += '"' + col_data[k] + '"' + delimiter
                        i += 1
                    else:
                        i += 4
                data += '\n'
        return title + "\n" + data
