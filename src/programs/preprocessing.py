import numpy as np


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
