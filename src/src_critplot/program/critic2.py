# -*- coding: utf-8 -*-
from typing import List
import copy
import os
import math
import numpy as np
import numpy.linalg

from core_gui_atomistic import helpers
from src_critplot.models.critical_point import CriticalPoint
from core_gui_atomistic.atom import Atom
from core_gui_atomistic.periodic_table import TPeriodTable
from src_critplot.models.atomic_model_cp import AtomicModelCP


def structure_from_cro_file(filename):
    """
    * Complete CP list, bcp and rcp connectivity table
    # (cp(end)+lvec connected to bcp/rcp)
    #cp  ncp   typ        position (cryst. coords.)            end1 (lvec)      end2 (lvec)
    1      1    n   0.41201875    0.19237040    0.32843532
    """
    #  print("structure_from_cro_file start")
    period_table = TPeriodTable()
    model = AtomicModelCP()
    if os.path.exists(filename) and filename.endswith("cro"):
        box_ang = helpers.from_file_property(filename, "Lattice parameters (ang):", 1, 'string').split()
        box_ang = np.array(helpers.list_str_to_float(box_ang))
        box_deg = helpers.from_file_property(filename, "Lattice angles (degrees):", 1, 'string').split()
        box_deg = np.array(helpers.list_str_to_float(box_deg))

        lat_vectors = helpers.lat_vectors_from_params(box_ang[0], box_ang[1], box_ang[2],
                                                      math.radians(box_deg[0]),
                                                      math.radians(box_deg[1]),
                                                      math.radians(box_deg[2]))

        model.set_lat_vectors(lat_vectors[0], lat_vectors[1], lat_vectors[2])

        crit_points = get_critical_points_info(filename)

        f = open(filename)
        str1 = f.readline()
        while str1.find("* Complete CP list, bcp and rcp connectivity table") < 0:
            str1 = f.readline()

        f.readline()
        f.readline()
        str1 = f.readline()
        while len(str1) > 5:
            str1 = str1.replace("(", "")
            str1 = str1.replace(")", "")
            data = str1.split()
            title_number = data[0]
            if data[0] != data[1]:
                title_number += "(" + data[1] + ")"
            number = int(data[1]) - 1
            """
            #cp  ncp   typ        position (cryst. coords.)            end1 (lvec)      end2 (lvec)
            1      1    n   0.67394965    0.51875056    0.41176494  
            115    37   b   0.61261657    0.49723994    0.38722411  101  ( 0  0  0 )  76  ( 0  0  0 )
            347    108  r   0.51856128    0.51875294    0.44655967  359  ( 0  0  0 ) 359  ( 0  0  0 )
            362    116  c   0.51852061    0.51875000    0.63729618
            """
            xyz = np.array([data[3], data[4], data[5]], dtype=float)
            crit_info = crit_points[number]
            if crit_info[3] == "nucleus":
                let = crit_info[8].replace("_", "")
                cp_type = "(3,-3)"
                title = let + title_number
                new_cp = init_crit_point(crit_info, let, cp_type, title, xyz)
                new_atom = Atom([*xyz, let, period_table.get_charge_by_letter(let)])
                model.add_atom(new_atom)
                # new_cp.let = "nn"
                model.add_critical_point(new_cp)

            if crit_info[3] == "nnattr":
                let = "attr"
                cp_type = "(3,-3)"
                title = let + title_number
                new_cp = init_crit_point(crit_info, let, cp_type, title, xyz)
                model.add_critical_point(new_cp)

            if crit_info[3] == "bond":
                let = "xb"
                cp_type = "(3,-1)"
                title = let + title_number
                new_cp = init_crit_point(crit_info, let, cp_type, title, xyz)
                add_atoms_to_cp(data, model, new_cp)
                model.add_critical_point(new_cp)

            if crit_info[3] == "ring":
                let = "xr"
                cp_type = "(3,+1)"
                title = let + title_number
                new_cp = init_crit_point(crit_info, let, cp_type, title, xyz)
                add_atoms_to_cp(data, model, new_cp)
                model.add_critical_point(new_cp)

            if crit_info[3] == "cage":
                let = "xc"
                cp_type = "(3,+3)"
                title = let + title_number
                new_cp = init_crit_point(crit_info, let, cp_type, title, xyz)
                add_atoms_to_cp(data, model, new_cp)
                model.add_critical_point(new_cp)
            str1 = f.readline()

        f.readline()
        f.close()
        model.convert_from_direct_to_cart()

        for cp in model.cps:
            ind1 = cp.get_property("atom1")
            ind2 = cp.get_property("atom2")
            if (ind1 is not None) and (ind2 is not None):
                trans1 = np.array(cp.get_property("atom1_translation"))
                trans2 = np.array(cp.get_property("atom2_translation"))
                p1 = CriticalPoint([model.cps[ind1 - 1].xyz + trans1, "xz", "bp"])
                p2 = CriticalPoint([cp.xyz, "xz", "bp"])
                p3 = CriticalPoint([model.cps[ind2 - 1].xyz + trans2, "xz", "bp"])
                if (np.linalg.norm(trans1) > 0) and (ind1 < model.n_atoms()):
                    atom = copy.deepcopy(model.atoms[ind1 - 1])
                    atom.xyz += trans1
                    atom.tag = "translated"
                    model.add_atom(atom, min_dist=-0.01)
                if (np.linalg.norm(trans2) > 0) and (ind2 < model.n_atoms()):
                    atom = copy.deepcopy(model.atoms[ind2 - 1])
                    atom.xyz += trans2
                    atom.tag = "translated"
                    model.add_atom(atom, min_dist=-0.01)
                model.add_bond_path_point([p2, p1])
                model.add_bond_path_point([p2, p3])
                atom_to_atom = model.cps[ind1 - 1].let + str(ind1) + "-" + model.cps[ind2 - 1].let + str(ind2)
                cp.set_property("atom_to_atom", atom_to_atom)
                bond_len = np.linalg.norm(model.cps[ind2 - 1].xyz + trans2 - model.cps[ind1 - 1].xyz - trans1)
                cp.set_property("cp_bp_len", bond_len)
        model.bond_path_points_optimize()
    return [model]


def add_atoms_to_cp(data, model, new_cp):
    if len(data) > 6:
        atom1 = int(data[6])
        atom2 = int(data[10])
        if (atom1 <= len(model.atoms)) and (atom2 <= len(model.atoms)):
            new_cp.set_property("atom1", atom1)
            new_cp.set_property("atom2", atom2)
            translation1 = int(data[7]) * model.lat_vector1 + int(data[8]) * model.lat_vector2 + \
                           int(data[9]) * model.lat_vector3
            translation2 = int(data[11]) * model.lat_vector1 + int(data[12]) * model.lat_vector2 + \
                           int(data[13]) * model.lat_vector3
            new_cp.set_property("atom1_translation", translation1)
            new_cp.set_property("atom2_translation", translation2)
        else:
            print("strange critical point: ", new_cp.to_string())


def get_critical_points_info(filename: str):
    f = open(filename)
    str1 = f.readline()
    while str1.find("* Critical point list, final report (non-equivalent cps)") < 0:
        str1 = f.readline()
    f.readline()
    f.readline()
    f.readline()
    str1 = f.readline()
    crit_points = []

    while len(str1) > 5:
        """
        ncp   pg  type   CPname         position (cryst. coords.)       mult  name            f             |grad|           lap
        1    C1  (3,-3) nucleus   0.41201875   0.19237040   0.32843532    1     H_      4.02547716E+00  0.00000000E+00 -1.04730465E+02
        """
        data1 = str1.split(")")
        data2 = data1[0].split("(")
        data = data1[1].split()
        row = data2[0].split()
        row.append(data2[1])
        row += data
        crit_points.append(row)
        str1 = f.readline()

    while (str1.find("Additional properties at the critical points") < 0) and (len(str1) > 0):
        str1 = f.readline()
    str1 = f.readline()
    point = 0

    while str1.find("+ Critical point no.") >= 0:
        text = ""
        str1 = f.readline()
        while not str1.startswith("+ "):
            if len(str1) > 0:
                text += str1
            str1 = f.readline()
        crit_points[point].append(text)
        point += 1
    f.close()
    return crit_points


def init_crit_point(crit_info, let, cp_type, title, xyz):
    new_cp = CriticalPoint([xyz, let, cp_type])
    new_cp.set_property("title", title)
    new_cp.set_property("rho", crit_info[9])
    new_cp.set_property("grad", crit_info[10])
    new_cp.set_property("lap", crit_info[11])
    new_cp.set_property("text", crit_info[12])
    return new_cp


def parse_bondpaths(filename: str, model: AtomicModelCP) -> List[AtomicModelCP]:
    """Import bond paths from *.xyz file.
    filename - name of file
    model - AtomicModelCP to add bondpaths from xyz file
    """
    molecules = []
    if os.path.exists(filename):
        f = open(filename)
        number_of_atoms = int(math.fabs(int(f.readline())))
        new_model = AtomicModelCP.atoms_from_xyz_structure(number_of_atoms, f, is_allow_charge_incorrect=True)
        model.delete_all_bond_paths()

        xz_points = []

        for atom in new_model.atoms:
            if atom.let.lower() == "xz":
                xz_points.append(atom)

        points = []

        for i in range(0, len(xz_points)):
            if len(points) == 0:
                points.append(xz_points[i])
            else:
                d = math.dist(points[-1].xyz, xz_points[i].xyz)
                if d < 0.09:
                    points.append(xz_points[i])
                else:
                    model.add_bond_path_point(points)
                    points = [xz_points[i]]

        if len(points) > 0:
            model.add_bond_path_point(points)

        model.bond_path_points_optimize()
        new_model = model
    molecules.append(new_model)
    return molecules


def model_to_critic_xyz_file(model, cps):
    """Returns data for *.xyz file with CP and BCP.
    model
    cps
    """
    text = ""

    n_atoms = model.n_atoms()
    for i in range(0, n_atoms):
        text += model.atoms[i].to_string() + "\n"

    n_cp = len(cps)
    for cp in cps:
        text += cp.to_string() + "\n"

    n_bcp = 0
    for cp in cps:
        bond1 = cp.bonds.get("bond1")
        bond2 = cp.bonds.get("bond2")

        for i in range(0, len(bond1)):
            n_bcp += 1
            text += bond1[i].to_string() + "\n"

        for i in range(0, len(bond2)):
            n_bcp += 1
            text += bond2[i].to_string() + "\n"

    header = "   " + str(n_atoms + n_cp + n_bcp) + "\n\n"
    return header + text


def create_critic2_xyz_file(bcp, bcp_seleсted, is_with_selected, model):
    if is_with_selected:
        text = model_to_critic_xyz_file(model, bcp_seleсted)
    else:
        for b in bcp_seleсted:
            for cp in bcp:
                if cp.to_string() == b.to_string():
                    bcp.remove(cp)
        text = model_to_critic_xyz_file(model, bcp)
    return text


def create_cri_file(cp_list, extra_points, is_form_bp, model, text_prop):
    sys_coord = np.array([model.lat_vector1, model.lat_vector2, model.lat_vector3])
    obr = np.linalg.inv(sys_coord).transpose()
    text = ""
    te = ""
    lines = ""
    textl = "crystal model.BADER.cube\n"
    textl += "WRITE model.xyz\n"
    textl += "load model.BADER.cube\n"
    textl += "load model.VT.DN.cube\n"
    textl += "load model.VT.UP.cube\n"
    textl += 'LOAD AS "-$2-$3"\n'
    textl += 'LOAD AS LAP 1\n'
    textl += "REFERENCE 1\n"

    for ind in cp_list:
        cp = model.cps[ind]
        text += "Bond Critical Point: " + str(ind) + "  :  "
        ind1 = cp.get_property("atom1")
        ind2 = cp.get_property("atom2")
        atom1 = model.atoms[ind1].let + str(ind1 + 1)
        atom2 = model.atoms[ind2].let + str(ind2 + 1)
        title = atom1 + "-" + atom2
        text += title + "\n"

        if is_form_bp:
            """ bond path """
            bond1 = cp.get_property("bond1")
            bond2 = cp.get_property("bond2")

            path_low = []
            for i in range(0, len(bond1)):
                index = len(bond1) - i
                coord = np.array([bond1[index - 1].x, bond1[index - 1].y, bond1[index - 1].z])
                res = obr.dot(coord)
                path_low.append(np.array([float(res[0]), float(res[1]), float(res[2])]))

            from_to = "{0:14.10} {1:14.10} {2:14.10} {3:14.10} {4:14.10} {5:14.10} ".format(path_low[0][0],
                                                                                            path_low[0][1],
                                                                                            path_low[0][2],
                                                                                            path_low[-1][0],
                                                                                            path_low[-1][1],
                                                                                            path_low[-1][2])

            lines += "# " + title + "\n"
            lines += "REFERENCE 1\n"
            lines += "LINE " + from_to + " 100 FILE ./lines/lines-" + title + "-00-charge.txt\n"
            lines += "REFERENCE 4\n"
            lines += "LINE " + from_to + " 100 FILE ./lines/lines-" + title + "-00-totpot.txt\n"
            lines += "REFERENCE 5\n"
            lines += "LINE " + from_to + " 100 FILE ./lines/lines-" + title + "-00-elpot.txt\n"

            first = "{0:14.10} {1:14.10} {2:14.10} ".format(path_low[-1][0], path_low[-1][1], path_low[-1][2])

            for i in range(1, len(bond2)):
                coord = np.array([bond2[i].x, bond2[i].y, bond2[i].z])
                res = obr.dot(coord)
                path_low.append(np.array([float(res[0]), float(res[1]), float(res[2])]))

            last = "{0:14.10} {1:14.10} {2:14.10} ".format(path_low[-1][0], path_low[-1][1], path_low[-1][2])
            lines += "REFERENCE 1\n"
            lines += "LINE " + first + last + " 100 FILE ./lines/lines-" + title + "-01-charge.txt\n"
            lines += "REFERENCE 4\n"
            lines += "LINE " + first + last + " 100 FILE ./lines/lines-" + title + "-01-totpot.txt\n"
            lines += "REFERENCE 5\n"
            lines += "LINE " + first + last + " 100 FILE ./lines/lines-" + title + "-01-elpot.txt\n"

            path_fine = [path_low[0]]
            for i in range(1, len(path_low)):
                dv = (path_low[i] - path_low[i - 1]) / extra_points
                for j in range(0, extra_points):
                    path_fine.append(path_fine[-1] + dv)

            for i in range(0, len(path_fine)):
                text += "{0:20.16f} {1:20.16f} {2:20.16f}\n".format(path_fine[i][0], path_fine[i][1], path_fine[i][2])
                te += "{0:20.16f} {1:20.16f} {2:20.16f}\n".format(path_fine[i][0], path_fine[i][1], path_fine[i][2])
        else:
            """ critical points only """
            res = obr.dot(np.array([cp.x, cp.y, cp.z]))
            text += "{0:20.16f} {1:20.16f} {2:20.16f}\n".format(float(res[0]), float(res[1]), float(res[2]))
            te += "{0:20.16f} {1:20.16f} {2:20.16f}\n".format(float(res[0]), float(res[1]), float(res[2]))
    textl += "# bond path information\n" + text_prop
    textl += 'POINTPROP elpot "$4"\n'
    textl += 'POINTPROP lapl "$5"\n'
    textl += "POINT ./POINTS.txt\n"
    lines += "UNLOAD ALL\nEND"
    return textl, lines, te, text
