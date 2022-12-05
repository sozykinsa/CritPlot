# -*- coding: utf-8 -*-

from core_gui_atomistic.atom import Atom


class AtomCp(Atom):
    """The atom class with bond paths."""

    def __init__(self, at_data):
        """Constructor"""
        super(AtomCp, self).__init__(at_data)
        self.bonds = {}
