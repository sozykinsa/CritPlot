# -*- coding: utf-8 -*-
import sys
import os
import math
import numpy as np
from core_gui_atomistic import helpers
from models.atomic_model_cp import AtomicModelCP
from core_gui_atomistic.periodic_table import TPeriodTable
sys.path.append('.')


def parse_cp_properties(filename, model):
    box_bohr, box_ang, box_deg, cps = check_cro_file(filename)
    al = math.radians(box_deg[0])
    be = math.radians(box_deg[1])
    ga = math.radians(box_deg[2])
    lat_vect_1, lat_vect_2, lat_vect_3 = helpers.lat_vectors_from_params(box_ang[0], box_ang[1], box_ang[2],
                                                                         al, be, ga)
    model.set_lat_vectors(lat_vect_1, lat_vect_2, lat_vect_3)
    k = 0
    for cp in model.bcp:
        for cp1 in cps:
            d1 = (cp.x - cp1[1]) ** 2
            d2 = (cp.y - cp1[2]) ** 2
            d3 = (cp.z - cp1[3]) ** 2
            if math.sqrt(d1 + d2 + d3) < 1e-5:
                cp.set_property("title", "b" + cp1[8])
                cp.set_property("field", cp1[4])
                cp.set_property("grad", cp1[5])
                cp.set_property("lap", cp1[6])
                cp.set_property("text", cp1[7])
                k += 1


def atoms_from_xyz(filename):
    """Import from *.xyz file."""
    molecules = []
    if os.path.exists(filename):
            f = open(filename)
            number_of_atoms = int(math.fabs(int(f.readline())))
            new_model = AtomicModelCP.atoms_from_xyz_structure(number_of_atoms, f)
            critic_data = {"xn", "xr", "xb", "xc", "xz"}
            fl = False
            for atom in new_model.atoms:
                if atom.let.lower() in critic_data:
                    fl = True
            if fl:
                    new_model2 = AtomicModelCP()
                    xz_points = []

                    for atom in new_model.atoms:
                        if atom.let.lower() not in critic_data:
                            new_model2.add_atom(atom)
                        if atom.let.lower() == "xb":
                            new_model2.add_critical_point_bond(atom)
                        if atom.let.lower() == "xz":
                            xz_points.append(atom)

                    points = []

                    for i in range(0, len(xz_points)):
                        if len(points) == 0:
                            points.append(xz_points[i])
                        else:
                            px = points[-1].x
                            py = points[-1].y
                            pz = points[-1].z

                            nx = xz_points[i].x
                            ny = xz_points[i].y
                            nz = xz_points[i].z

                            d = math.sqrt((px - nx) * (px - nx) + (py - ny) * (py - ny) + (pz - nz) * (pz - nz))

                            if d < 0.09:
                                points.append(xz_points[i])
                            else:
                                new_model2.add_bond_path_point(points)
                                points = [xz_points[i]]

                    if len(points) > 0:
                        new_model2.add_bond_path_point(points)

                    new_model2.bond_path_points_optimize()
                    new_model = new_model2
    molecules.append(new_model)
    return molecules


def check_cro_file(filename):
    if os.path.exists(filename) and filename.endswith("cro"):
        box_bohr = helpers.from_file_property(filename, "Lattice parameters (bohr):", 1, 'string').split()
        box_bohr = np.array(helpers.list_str_to_float(box_bohr))
        box_ang = helpers.from_file_property(filename, "Lattice parameters (ang):", 1, 'string').split()
        box_ang = np.array(helpers.list_str_to_float(box_ang))
        box_deg = helpers.from_file_property(filename, "Lattice angles (degrees):", 1, 'string').split()
        box_deg = np.array(helpers.list_str_to_float(box_deg))

        filename = open(filename)
        str1 = filename.readline()
        while (str1.find("Critical point list, final report (non-equivalent cps") < 0) and (len(str1) > 0):
            str1 = filename.readline()
        filename.readline()
        filename.readline()
        filename.readline()

        cps = []
        str1 = filename.readline()

        while len(str1) > 3:
            str1 = str1.split(')')[1].split()
            x = float(str1[1]) * box_ang[0]
            y = float(str1[2]) * box_ang[1]
            z = float(str1[3]) * box_ang[2]

            line = [str1[0], x, y, z, str1[6], str1[7], str1[8], "", ""]
            cps.append(line)
            str1 = filename.readline()

        while (str1.find("Additional properties at the critical points") < 0) and (len(str1) > 0):
            str1 = filename.readline()
        str1 = filename.readline()
        point = 0
        while str1.find("+ Critical point no.") >= 0:
            text = ""
            title = str1.split("+ Critical point no.")
            str1 = filename.readline()
            while not str1.startswith("+ "):
                if len(str1) > 0:
                    text += str1
                str1 = filename.readline()
            cps[point][7] = text
            cps[point][8] = cps[point][0][0] + helpers.spacedel(title[1])
            point += 1

        filename.close()
        return box_bohr, box_ang, box_deg, cps
    else:
        return "", "", "", []


def create_csv_file_cp(cp_list, model, delimiter: str = ";"):
    title = ""
    data = ""
    for ind in cp_list:
        title = ""
        cp = model.bcp[ind]
        title += "BCP" + delimiter + "atoms" + delimiter + "dist" + delimiter
        data += str(ind) + delimiter
        ind1 = cp.get_property("atom1")
        ind2 = cp.get_property("atom2")
        atom1 = model.atoms[ind1].let + str(ind1 + 1)
        atom2 = model.atoms[ind2].let + str(ind2 + 1)
        data += atom1 + "-" + atom2 + delimiter

        dist_line = round(model.atom_atom_distance(ind1, ind2), 4)
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


def model_to_critic_xyz_file(model, cps):
    """Returns data for *.xyz file with CP and BCP."""
    text = ""

    n_atoms = model.n_atoms()
    for i in range(0, n_atoms):
        text += model.atoms[i].to_string() + "\n"

    n_cp = len(cps)
    for cp in cps:
        text += cp.to_string() + "\n"

    n_bcp = 0
    for cp in cps:
        bond1 = cp.get_property("bond1")
        bond2 = cp.get_property("bond2")

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
        cp = model.bcp[ind]
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
