# -*- coding: utf-8 -*-

from src_core_atomistic.atom import Atom


class CriticalPoint(Atom):
    """The atom class with bond paths."""

    def __init__(self, at_data):
        """Constructor"""
        self.xyz = at_data[0]
        self.let = at_data[1]
        self.cp_type: str = at_data[2]
        self.charge = 0
        self.is_visible: bool = True
        self.selected: bool = False
        self.active: bool = False
        self.fragment1: bool = False
        self.properties = {}
        self.visible_property = ""
        self.tag = ""
        self.bonds = {}
