# This file is a part of CritPlot programm.

from typing import Callable

import OpenGL.GL as gl
import OpenGL.GLU as glu
import numpy as np
from core_gui_atomistic.guiopenglbase import GuiOpenGLBase
from core_gui_atomistic.helpers import is_number
from src_critplot.models.atomic_model_cp import AtomicModelCP


class GuiOpenGLCP(GuiOpenGLBase):
    def __init__(self, parent=None):
        super().__init__(parent)

        # opengl lists
        self.list_for_bcp: int = 8
        self.list_for_bondpath: int = 9

        self.is_bcp_available: bool = False
        self.is_bond_path_available: bool = False
        self.is_show_bcp: bool = True
        self.is_show_bond_path: bool = True
        self.selected_cp = -1
        self.is_bcp_property_visible: bool = False
        self.bcp_property: str = ""
        self.NLists = 10

        self.selected_cp_callback: Callable = None

        self.main_model = AtomicModelCP()

    def init_params(self, the_object) -> None:
        super().init_params(the_object)
        self.selected_cp = the_object.selected_cp

    def add_all_elements(self):
        super().add_all_elements()
        self.add_critical_points()
        self.add_bond_path()

    def set_form_elements(self, check_atom_selection=None, orientation_model_changed: Callable = None,
                          selected_atom_position: Callable = None, selected_atom_changed: Callable = None,
                          selected_cp_changed: Callable = None, quality=1):
        """Set pointers for Form update.
            Args:
                check_atom_selection: ...
                orientation_model_changed: pointer to MainForm.orientation_model_changed();
                selected_atom_position: pointer to MainForm.selected_atom_position();
                selected_atom_changed: pointer to MainForm.selected_atom_changed();
                selected_cp_changed: pointer to MainForm.selected_cp_changed().
                quality: ... .
        """
        super().set_form_elements(check_atom_selection, orientation_model_changed, selected_atom_position,
                                  selected_atom_changed, quality)
        self.selected_cp_callback = selected_cp_changed

    def set_property_show_bcp(self, show_bcp):
        self.is_show_bcp = show_bcp

    def set_property_bond_path(self, bond_path):
        self.is_show_bond_path = bond_path

    def show_bcp_property(self, prop: str = ""):
        if prop == "":
            self.is_bcp_property_visible = False
        else:
            self.bcp_property = prop
            for bcp in self.main_model.bcp:
                value = bcp.get_property(prop)
                if is_number(value):
                    value = float(value)
                    value = str(round(value, self.property_precision))
                bcp.visible_property = value
            self.is_bcp_property_visible = True
        self.update()

    def init_params(self, ogl_model) -> None:
        super().init_params(ogl_model)
        self.selected_cp = ogl_model.selected_cp
        self.is_bcp_property_visible = ogl_model.is_bcp_property_visible
        self.bcp_property = ogl_model.bcp_property
        self.is_bcp_available = ogl_model.is_bcp_available
        self.is_bond_path_available = ogl_model.is_bond_path_available
        self.is_show_bcp = ogl_model.is_show_bcp
        self.is_show_bond_path = ogl_model.is_show_bond_path

    def copy_state(self, ogl_model):
        super().copy_state(ogl_model)
        self.add_critical_points()
        self.add_bond_path()
        self.update()

    def set_structure_parameters(self, atoms_colors, is_view_atoms, is_view_atom_numbers, is_view_box, box_color,
                                 is_view_bonds, bonds_color, bond_width, bonds_by_atoms, is_view_axes, axes_color,
                                 contour_width, is_bcp_property_visible):
        super().set_structure_parameters(atoms_colors, is_view_atoms, is_view_atom_numbers, is_view_box,
                                         box_color, is_view_bonds, bonds_color, bond_width, bonds_by_atoms,
                                         is_view_axes, axes_color, contour_width)
        self.is_bcp_property_visible = is_bcp_property_visible
        self.add_critical_points()
        self.add_bond_path()
        self.update()

    def selected_atom_properties_to_form(self):
        super().selected_atom_properties_to_form()
        self.selected_cp_callback(self.selected_cp)

    def add_critical_points(self):
        gl.glNewList(self.object + 8, gl.GL_COMPILE)
        color = (1, 0, 0)
        self.add_cp(color, self.main_model.bcp)

        color = (1, 1, 0)
        self.add_cp(color, self.main_model.ccp)

        color = (1, 1, 1)
        self.add_cp(color, self.main_model.rcp)

        gl.glEndList()
        self.is_bcp_available = True
        self.update()

    def add_cp(self, color, cps):
        for at in cps:
            gl.glPushMatrix()
            gl.glTranslatef(*(self.scale_factor * at.xyz))
            gl.glColor3f(*color)
            mult = self.scale_factor
            if at.is_selected():
                gl.glColor3f(0, 0, 1)
                mult *= 1.3
            glu.gluSphere(glu.gluNewQuadric(), 0.15 * mult, self.quality * 70, self.quality * 70)
            gl.glPopMatrix()

    def add_bond_path(self):
        gl.glNewList(self.object + 9, gl.GL_COMPILE)

        for cp in self.main_model.bcp:
            self.add_critical_path(cp.get_property("bond1opt"))
            self.add_critical_path(cp.get_property("bond2opt"))

        gl.glEndList()
        self.is_bond_path_available = True
        self.update()

    def add_critical_path(self, bond):
        if not bond:
            return

        if np.linalg.norm(bond[0].xyz - bond[-1].xyz) > 4.0:
            return

        gl.glColor3f(0, 1, 0)
        for i in range(1, len(bond)):
            self.add_bond(self.scale_factor * bond[i - 1].xyz,
                          self.scale_factor * bond[i].xyz,
                          self.scale_factor * 0.03)

    def paintGL(self):
        self.makeCurrent()
        self.configure_gl()
        try:
            self.prepere_scene()
            self.light_prepare()
            if self.active:
                self.prepare_orientation()
                if self.is_view_atoms:
                    gl.glCallList(self.object + self.list_for_atoms)  # atoms

                if self.is_bcp_available and self.is_show_bcp:
                    gl.glCallList(self.object + 8)  # BCP

                if self.can_atom_search:
                    self.get_atom_on_screen()

                if self.is_view_bonds and (len(self.main_model.bonds) > 0):
                    gl.glCallList(self.object + 2)  # find_bonds_exact

                if self.is_view_box:
                    gl.glCallList(self.object + 3)  # lattice_parameters_abc_angles

                if self.is_view_axes:
                    gl.glCallList(self.object + 7)  # Axes

                if self.is_bond_path_available and self.is_show_bond_path:
                    gl.glCallList(self.object + 9)  # Bondpath

                text_to_render = []
                if self.is_atomic_numbers_visible:
                    for i in range(0, len(self.main_model.atoms)):
                        at = self.main_model.atoms[i]
                        text_to_render.append([self.scale_factor * at.x, self.scale_factor * at.y,
                                               self.scale_factor * at.z, at.let + str(i + 1)])

                if self.is_bcp_property_visible:
                    for i in range(0, len(self.main_model.bcp)):
                        at = self.main_model.bcp[i]
                        text_to_render.append([self.scale_factor * at.x, self.scale_factor * at.y,
                                               self.scale_factor * at.z, at.visible_property])

                if self.is_atomic_numbers_visible or self.is_bcp_property_visible:
                    self.render_text(text_to_render)
        except Exception as exc:
            print(exc)
            pass

    def get_atom_on_screen(self):
        point = self.get_point_in_3d(self.x_scene, self.y_scene)
        old_selected = self.selected_atom
        old_selected_cp = self.selected_cp
        need_for_update = False

        atom_ind, atom_min_r = self.nearest_point(self.scale_factor, self.main_model.atoms, point)
        bcp_ind, bcp_min_r = self.nearest_point(self.scale_factor, self.main_model.bcp, point)
        ccp_ind, ccp_min_r = self.nearest_point(self.scale_factor, self.main_model.ccp, point)
        rcp_ind, rcp_min_r = self.nearest_point(self.scale_factor, self.main_model.rcp, point)

        if (bcp_min_r < 1.4) and (bcp_min_r <= atom_min_r):
            if self.selected_cp == bcp_ind:
                if self.selected_cp > 0:
                    self.main_model.bcp[self.selected_cp].set_selected(False)
                self.selected_cp = -1
            else:
                if self.selected_cp > 0:
                    self.main_model.bcp[self.selected_cp].set_selected(False)
                self.selected_cp = bcp_ind
                if self.selected_atom > 0:
                    self.main_model.atoms[self.selected_atom].set_selected(False)
                self.selected_atom = -1
                if self.selected_cp > 0:
                    self.main_model.bcp[self.selected_cp].set_selected(True)

        if atom_min_r <= bcp_min_r:
            self.update_selected_atom(atom_ind, atom_min_r)

        self.can_atom_search = False
        selected_atom_was_modified = old_selected != self.selected_atom
        selected_cp_was_modified = old_selected_cp != self.selected_cp
        if selected_atom_was_modified:
            if self.selected_cp > 0:
                self.main_model.bcp[self.selected_cp].set_selected(False)
                self.selected_cp = -1
        if selected_atom_was_modified or selected_cp_was_modified:
            self.add_atoms()
            self.add_critical_points()
            self.add_bonds()
            self.selected_atom_changed()
            self.update()