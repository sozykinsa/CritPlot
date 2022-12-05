# -*- coding: utf-8 -*-
import os
import numpy as np
from core_gui_atomistic import helpers
from core_gui_atomistic.periodic_table import TPeriodTable
from src_critplot.models.atom_cp import AtomCp as Atom
from src_critplot.models.atomic_model_cp import AtomicModelCP


def number_of_atoms_from_outp(filename):
    """Number of atoms in TOPOND output."""
    n_atoms = 0
    if os.path.exists(filename):
        result = helpers.from_file_property(filename, 'N. OF ATOMS PER CELL', prop_type='string')
        n_atoms = int(result.split()[0])
    return n_atoms


def get_cell(filename):
    """Lattice vectors."""
    if os.path.exists(filename):
        k = -1
        out = 0
        with open(filename, "r") as file1:
            lat_vectors = 500 * np.eye(3)
            for line in file1:
                if k >= 0:
                    data = line.split()
                    if len(data) == 3:
                        lat_vectors[3-k] = np.array(data)
                    k -= 1
                if line.find("DIRECT LATTICE VECTOR COMPONENTS (ANGSTROM)") >= 0:
                    k = 3
                    out = 1
                if (k == 0) and (out == 1):
                    return lat_vectors
    return None


def get_atoms(filename):
    """Atoms positions."""
    model = AtomicModelCP()
    period_table = TPeriodTable()
    if os.path.exists(filename):
        number_of_atoms = number_of_atoms_from_outp(filename)
        file1 = open(filename)
        row = file1.readline()
        while row and (row.find("ATOM N.AT.  SHELL    X(A)      Y(A)      Z(A)      EXAD       N.ELECT.") < 0):
            row = file1.readline()
        if row:
            file1.readline()
            for n in range(number_of_atoms):
                data = file1.readline().split()
                charge = period_table.get_charge_by_letter(data[2])
                model.add_atom(Atom([float(data[4]), float(data[5]), float(data[6]), data[2], charge]))
    return model


def parse_cp_data(filename: str, model: AtomicModelCP):
    """Get critical points data."""
    if os.path.exists(filename):
        file1 = open(filename)
        row = file1.readline()
        while row:
            while row and (row.find("CP N.") < 0):
                row = file1.readline()
            title = row.split("CP N.")[1].split()[0]
            row1 = row.split()

            file1.readline()
            file1.readline()
            data1 = file1.readline().split(":")
            data = helpers.spacedel(data1[1])

            if data == "(3,-1)":
                """
                CP N.      8  NON-EQUIV. ATOM      3 C  ---    199 H (     1      1      0 )  DISTANCE(ANG)= 1.0859
                *********

                CP TYPE                        :  (3,-1)
                COORD(AU)  (X  Y  Z)           :  1.8343E+00  7.5829E+00  2.5685E+01
                COORD FRACT. CONV. CELL        :  3.5525E-02  1.4686E-01  4.9744E-01
                PROPERTIES (RHO,GRHO,LAP)      :  2.7486E-01  1.2022E-15 -8.2652E-01
                KINETIC ENERGY DENSITIES (G,K) :  3.7054E-02  2.4369E-01
                VIRIAL DENSITY                 : -2.8074E-01
                ELF(PAA)                       :  9.8781E-01
                """
                text = "Type : (3,-1)\n"
                title = "b" + title
                let = "xb"
                cp, row = parse_cp_point(file1, let, text, title)
                add_associated_atoms(cp, model, row1)
                model.add_critical_point(cp)
            elif data == "(3,-3)":
                """
                CP N.      9
                **********

                ATTRACTOR CP TYPE              :  (3,-3)
                COORD(AU)  (X  Y  Z)           :  2.2737E+00  7.4470E+00  2.6232E+01
                COORD FRACT. CONV. CELL        :  4.4035E-02  1.4423E-01 -4.9195E-01
                PROPERTIES (RHO,GRHO,LAP)      :  4.3196E-01  3.5396E-14 -2.1931E+01
                TRAJECTORY LENGTH(ANG)         :  3.7842E-01
                INTEGRATION STEPS              :      21
                """
                # print("nucleus (3,-3) ----->")
                text = "Type : (3,-3)\n"
                title = "A" + title
                let = "A"
                cp, row = parse_cp_point(file1, let, text, title)
                model.add_critical_point(cp)
                # print(data)
            elif data == "(3,+1)":
                """ CP TYPE                        :  (3,+1) """
                # print("ring (3,+1) ----->")
                """
                CP N.      7  NON-EQUIV. ATOM      2 O  ---     31 O (     0      0      0 )  DISTANCE(ANG)= 2.3736
                *********

                CP TYPE                        :  (3,+1)
                COORD(AU)  (X  Y  Z)           : -1.7694E+00  1.7694E+00 -2.4047E+01
                COORD FRACT. CONV. CELL        : -3.4268E-02  3.4268E-02 -4.6573E-01
                PROPERTIES (RHO,GRHO,LAP)      :  2.7570E-02  1.4387E-16  1.3993E-01
                KINETIC ENERGY DENSITIES (G,K) :  3.4069E-02 -9.1500E-04
                VIRIAL DENSITY                 : -3.3154E-02
                ELF(PAA)                       :  4.3030E-02
                """

                text = "Type : (3,+1)\n"
                title = "r" + title
                let = "xr"
                cp, row = parse_cp_point(file1, let, text, title)
                add_associated_atoms(cp, model, row1)
                model.add_critical_point(cp)
            elif data == "(3,+3)":
                """ CP TYPE                        :  (3,+3) """
                # print("cage (3,+3) ----->")
                text = "Type : (3,+3)\n"
                title = "c" + title
                let = "xc"
                cp, row = parse_cp_point(file1, let, text, title)
                add_associated_atoms(cp, model, row1)
                model.add_critical_point(cp)
                # print(data)
            else:
                print("else")
                print(row)
                row = file1.readline()
            while (row.find("NUMBER OF UNIQUE CRI. POINT FOUND") < 0) and (row.find("CP N.") < 0):
                row = file1.readline()
            if row.find("NUMBER OF UNIQUE CRI. POINT FOUND") > 0:
                """ correction """
                n_cp = int(helpers.spacedel(row.split("NUMBER OF UNIQUE CRI. POINT FOUND:")[1]))
                # print("n_cp: ", n_cp)
                row = file1.readline()
                while len(row) < 10:
                    row = file1.readline()
                row = file1.readline()
                while len(row) < 10:
                    row = file1.readline()
                """CP N.  X(AU)   Y(AU)   Z(AU)   TYPE      RHO    LAPL   L1(V1) L2(V2) L3(V3) ELLIP"""
                while row.find("*************************************") < 0:
                    row = file1.readline()
                for i in range(n_cp):
                    file1.readline()
                    row1 = file1.readline()
                    if len(row1) < 4:
                        row1 = file1.readline().split()
                        row2 = file1.readline().split()
                        # print("1---> ", row1)
                        # print("2---> ", row2)
                        b_len = round((float(row1[6]) + float(row2[6])) * 0.52917720859, 4)

                        cp = model.cps[i]
                        atom1_old = cp.get_property("atom1")
                        atom2_old = cp.get_property("atom2")
                        atom1 = int(row1[2])
                        atom2 = int(row2[2])
                        correction = "-"
                        text_new = cp.get_property("text")
                        f1 = (atom1_old == atom1 - 1) and (atom2_old == atom2 - 1)
                        if not (f1 or (atom1_old == atom2 - 1) and (atom2_old == atom1 - 1)):
                            cp.set_property("atom1", atom1 - 1)
                            cp.set_property("atom2", atom2 - 1)
                            let1 = model.atoms[atom1_old].let
                            let2 = model.atoms[atom2_old].let
                            dist_line = round(model.atom_atom_distance(atom1_old, atom2_old), 4)
                            let1n = model.atoms[atom1 - 1].let
                            let2n = model.atoms[atom2 - 1].let
                            correction = "(" + let1 + str(atom1_old + 1) + "-" + let2 + str(atom2_old + 1) + ")r=" + \
                                         str(dist_line) + "->(" + let1n + str(atom1) + "-" + let2n + str(atom2) + ")" +\
                                         "r=" + str(b_len)
                        text_new += "\ncorrections : " + correction
                        cp.set_property("text", text_new)
                        for j in range(3):
                            file1.readline()
                    else:
                        for j in range(4):
                            row = file1.readline()
                row = ""
        # print("cp -path test")
        for cp in model.cps:
            ind1 = cp.get_property("atom1")
            ind2 = cp.get_property("atom2")
            if (ind1 is not None) and (ind2 is not None):
                p1 = Atom([*model.atoms[ind1].xyz, "xz", 1])
                p2 = Atom([*cp.xyz, "xz", 1])
                p3 = Atom([*model.atoms[ind2].xyz, "xz", 1])
                model.add_bond_path_point([p2, p1])
                model.add_bond_path_point([p2, p3])
        model.bond_path_points_optimize()


def parse_cp_point(file1, let, text, title):
    row = file1.readline()
    text += helpers.spacedel(row) + "\n"
    """COORD(AU)  (X  Y  Z)           :  3.3821E+00  2.2164E+00 -3.7940E-16"""
    data = row.split(":")[1].split()
    x = float(data[0]) * 0.52917720859
    y = float(data[1]) * 0.52917720859
    z = float(data[2]) * 0.52917720859
    cp = Atom([x, y, z, let, 1])
    cp.set_property("title", title)
    text += helpers.spacedel(file1.readline()) + "\n"
    """COORD FRACT. CONV. CELL        :  2.5000E-01  2.5000E-01 -2.0484E-17"""
    row = helpers.spacedel(file1.readline()) + "\n"
    """PROPERTIES (RHO,GRHO,LAP)      :  3.1585E-03  2.9929E-18  1.3881E-02"""
    data = row.split(":")[1].split()
    cp.set_property("rho", data[0])
    text += "rho : " + data[0] + "\n"
    cp.set_property("grad", data[1])
    text += "grad : " + data[1] + "\n"
    cp.set_property("lap", data[2])
    text += "lap : " + data[2] + "\n"
    row = file1.readline()
    if (let == "xb") or (let == "xr"):
        """KINETIC ENERGY DENSITIES (G,K) :  2.4448E-03 -1.0254E-03"""
        text += "KINETIC ENERGY DENSITIES (G) : " + row.split()[5] + "\n"
        text += "KINETIC ENERGY DENSITIES (K) : " + row.split()[6] + "\n"
        row = helpers.spacedel(file1.readline())
        """VIRIAL DENSITY                 : -1.4194E-03"""
        text += row + "\n"
        row = helpers.spacedel(file1.readline()) + "\n"
        """ELF(PAA)                       :  6.3365E-03"""
        text += helpers.spacedel(row) + "\n"
        for i in range(8):
            file1.readline()
        row = file1.readline()
        """ELLIPTICITY                    :  1.4410E-01"""
        text += helpers.spacedel(row)
    cp.set_property("text", text)
    return cp, row


def add_associated_atoms(cp, model, row1):
    ind1, ind2 = int(row1[5]) - 1, int(row1[8]) - 1
    cp.set_property("atom1", ind1)
    cp.set_property("atom2", ind2)
    # print(row1)
    translation1 = np.zeros(3, dtype=float)
    translation2 = int(row1[11]) * model.lat_vector1 + int(row1[12]) * model.lat_vector2 + \
                   int(row1[13]) * model.lat_vector3
    cp.set_property("atom1_translation", translation1)
    cp.set_property("atom2_translation", translation2)


def atomic_data_from_output(filename):
    """import lattice and positions from TOPOND output."""
    model = AtomicModelCP()
    if os.path.exists(filename):
        lat_vectors = get_cell(filename)
        model = get_atoms(filename)
        model.set_lat_vectors(lat_vectors[0], lat_vectors[1], lat_vectors[2])
        parse_cp_data(filename, model)
    return [model]
