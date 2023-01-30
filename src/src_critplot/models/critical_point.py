# -*- coding: utf-8 -*-

from core_gui_atomistic.atom import Atom


class CriticalPoint(Atom):
    """The atom class with bond paths."""

    def __init__(self, at_data):
        """Constructor"""
        super().__init__(at_data)
        self.bonds = {}
