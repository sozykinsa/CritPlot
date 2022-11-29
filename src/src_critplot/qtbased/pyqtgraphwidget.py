# -*- coding: utf-8 -*-
# ------------------------------------------------------
# ------------------ PyqtGraphWidget -------------------
# ------------------------------------------------------
# https://www.pythonguis.com/tutorials/pyside-plotting-pyqtgraph/

import pyqtgraph as pg  # pip install pyqtgraph
import numpy as np
from typing import List
from core_gui_atomistic.pyqtgraphwidget_core import PyqtGraphWidgetCore


class PyqtGraphWidget(PyqtGraphWidgetCore):

    COLORS = [(0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255)]

    def __init__(self, parent=None):
        super().__init__(parent)

    def plot(self, x: List[List[float]], y: List[List[float]], labels: List[str],
             title: str, x_title: str, y_title: str, is_colored=True):
        self.title = title
        self.x_title = x_title
        self.y_title = y_title
        self.apply_styles()
        self.add_title()
        self.add_axes_titles()

        n_plots = len(y)
        for index in range(n_plots):
            pen = pg.mkPen(color=self.COLORS[index % len(self.COLORS) if is_colored else 0], width=self.line_width)
            self.graphWidget.plot(x[index % len(x)], y[index], name=labels[index % len(labels)],
                                  pen=pen, font=self.legend_font)

    def add_title(self):
        self.graphWidget.setTitle(self.title, color=self.font_color, size=str(self.font_size_title) + "pt")

    def add_scatter(self, xs, ys):
        scatter = pg.ScatterPlotItem(size=15, brush=pg.mkBrush(255, 255, 0, 190))
        spots = [{'pos': [xs[i], ys[i]], 'data': 1} for i in range(len(xs))]
        scatter.addPoints(spots)
        self.graphWidget.addItem(scatter)

    def add_axes_titles(self):
        # Add Axis Labels
        self.graphWidget.setLabel("left", self.y_title, **self.styles)
        self.graphWidget.setLabel("bottom", self.x_title, **self.styles)

    def add_line(self, _pos, _angle, _width, _style):  # pragma: no cover
        # _style: Qt.SolidLine, Qt.DashLine, Qt.DotLine, Qt.DashDotLine, Qt.DashDotDotLine
        pen = pg.mkPen(color=self.COLORS[0], width=_width, style=_style)
        line = pg.InfiniteLine(pos=_pos, angle=_angle, pen=pen)
        self.graphWidget.addItem(line)

    def add_histogram(self, values, num_bins, face_color, title, x_title, y_title) -> None:  # pragma: no cover
        self.apply_styles()
        y, x = np.histogram(values, bins=num_bins)
        curve = pg.PlotCurveItem(x, y, stepMode=True, fillLevel=0, brush=face_color)
        self.graphWidget.addItem(curve)
        self.x_title = x_title
        self.y_title = y_title
        self.title = title
        self.add_axes_titles()
        self.add_title()
