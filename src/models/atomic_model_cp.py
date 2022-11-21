# -*- coding: utf-8 -*-
# Python 3

import math
from copy import deepcopy
from numpy.linalg import inv
from core_gui_atomistic.atom import Atom
from core_gui_atomistic.atomic_model import AtomicModel


class AtomicModelCP(AtomicModel):
    def __init__(self, new_atoms: list = []):
        super().__init__(new_atoms)
        self.bcp = []

    def bond_path_points_optimize(self):
        i = 0

        while i < len(self.bcp):
            bond1 = self.bcp[i].get_property("bond1")
            bond2 = self.bcp[i].get_property("bond2")
            if (bond1 is None) or (bond2 is None):
                self.bcp.pop(i)
                i -= 1
            else:
                if (bond1[-1].x == bond2[-1].x) and (bond1[-1].y == bond2[-1].y) and (bond1[-1].z == bond2[-1].z):
                    self.bcp.pop(i)
                    i -= 1
            i += 1

        for cp in self.bcp:
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

        for point in self.bcp:
            point.x += l_x
            point.y += l_y
            point.z += l_z

            self.move_bond_path(l_x, l_y, l_z, point.get_property("bond1"))
            self.move_bond_path(l_x, l_y, l_z, point.get_property("bond2"))

            self.move_bond_path(l_x, l_y, l_z, point.get_property("bond1opt"))
            self.move_bond_path(l_x, l_y, l_z, point.get_property("bond2opt"))
        return self.atoms

    @staticmethod
    def move_bond_path(lx, ly, lz, bond):
        if bond:
            for bp in bond:
                bp.x += lx
                bp.y += ly
                bp.z += lz

    def go_to_positive_coordinates_translate(self):
        self.go_to_positive_array_translate(self.atoms)
        self.go_to_positive_array_translate(self.bcp)
        self.bond_path_opt_update()

    def bond_path_opt_update(self):
        for i in range(len(self.bcp)):
            self.bcp[i].properties.pop("bond1")
            self.bcp[i].properties.pop("bond2")
            self.bcp[i].properties.pop("bond1opt")
            self.bcp[i].properties.pop("bond2opt")

            ind1 = self.bcp[i].get_property("atom1")
            ind2 = self.bcp[i].get_property("atom2")
            p1 = Atom([*self.atoms[ind1].xyz, "xz", 1])
            p2 = Atom([*self.bcp[i].xyz, "xz", 1])
            p3 = Atom([*self.atoms[ind2].xyz, "xz", 1])

            self.add_bond_path_point([p2, p1])
            self.add_bond_path_point([p2, p3])
        self.bond_path_points_optimize()

    def n_bcp(self):
        return len(self.bcp)

    def move_atoms_to_cell(self):
        a_inv = inv(self.lat_vectors)
        self.move_object_to_cell(self.atoms, a_inv)
        self.move_object_to_cell(self.bcp, a_inv)

    def go_to_positive_coordinates(self):
        xm = self.minX()
        ym = self.minY()
        zm = self.minZ()
        self.go_to_positive_array(self.atoms, xm, ym, zm)
        self.go_to_positive_array(self.bcp, xm, ym, zm)

    def add_critical_point_bond(self, atom):
        new_at = deepcopy(atom)
        self.bcp.append(new_at)

    def add_bond_path_point(self, points):
        for cp in self.bcp:
            dx = math.pow(cp.x - points[0].x, 2)
            dy = math.pow(cp.y - points[0].y, 2)
            dz = math.pow(cp.z - points[0].z, 2)
            d = math.sqrt(dx + dy + dz)
            if d < 1e-4:
                if cp.get_property("bond1") is None:
                    cp.set_property("bond1", deepcopy(points))
                else:
                    cp.set_property("bond2", deepcopy(points))
                    ind1, ind2 = self.atoms_of_bond_path(self.atoms, self.bcp[self.bcp.index(cp)])
                    cp.set_property("atom1", ind1)
                    cp.set_property("atom2", ind2)

    @staticmethod
    def atoms_of_bond_path(atoms, cp):
        bond1 = cp.get_property("bond1")
        bond2 = cp.get_property("bond2")
        if bond1 is None or bond2 is None:
            return 0, 0
        cpx1 = bond1[-1].x
        cpy1 = bond1[-1].y
        cpz1 = bond1[-1].z
        cpx2 = bond2[-1].x
        cpy2 = bond2[-1].y
        cpz2 = bond2[-1].z
        x1 = atoms[0].x
        y1 = atoms[0].y
        z1 = atoms[0].z
        minr1 = math.sqrt((cpx1 - x1) * (cpx1 - x1) + (cpy1 - y1) * (cpy1 - y1) + (cpz1 - z1) * (cpz1 - z1))
        minr2 = math.sqrt((cpx2 - x1) * (cpx2 - x1) + (cpy2 - y1) * (cpy2 - y1) + (cpz2 - z1) * (cpz2 - z1))
        ind1 = 0
        ind2 = 0
        for i in range(0, len(atoms)):
            x1 = atoms[i].x
            y1 = atoms[i].y
            z1 = atoms[i].z

            dx1 = cpx1 - x1
            dx2 = cpx2 - x1

            dy1 = cpy1 - y1
            dy2 = cpy2 - y1

            dz1 = cpz1 - z1
            dz2 = cpz2 - z1

            d1 = math.sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1)
            d2 = math.sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2)

            if d1 < minr1:
                minr1 = d1
                ind1 = i

            if d2 < minr2:
                minr2 = d2
                ind2 = i
        return ind1, ind2
