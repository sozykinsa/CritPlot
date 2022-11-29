# -*- coding: utf-8 -*-
# ------------------------------------------------------
# ------------------ PyqtGraphWidget -------------------
# ------------------------------------------------------
# https://www.pythonguis.com/tutorials/pyside-plotting-pyqtgraph/

from PySide2.QtWidgets import QWidget, QVBoxLayout
from PySide2.QtGui import QFont
import pyqtgraph as pg  # pip install pyqtgraph
import numpy as np
from typing import List


class PyqtGraphWidgetCore(QWidget):

    COLORS = [(0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255)]

    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.graphWidget = pg.PlotWidget()
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.graphWidget)
        self.setLayout(vertical_layout)

        self.line_width = 2

        self.font_size_title = 20
        self.font_size_axes = 20
        self.font_size_legend = 20

        self.font_color = self.COLORS[0]
        self.font_pen = pg.mkPen(color=self.font_color)

        self.graphWidget.setBackground((255, 255, 255))
        self.legend_item = None

        self.styles = {"color": "#000", "font-size": str(self.font_size_axes) + "px"}

        self.title_font = QFont()
        self.legend_font = QFont()
        self.axes_font = QFont()

        self.title = ""
        self.x_title = ""
        self.y_title = ""

        self.apply_styles()

    def set_styles(self, font_size_t, font_size_a, font_size_l, line_width, font_color):
        self.font_size_title = font_size_t
        self.font_size_axes = font_size_a
        self.font_size_legend = font_size_l
        self.line_width = line_width
        self.font_color = font_color
        self.apply_styles()
