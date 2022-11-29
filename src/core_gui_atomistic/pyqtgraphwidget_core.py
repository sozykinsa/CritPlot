# -*- coding: utf-8 -*-
# ------------------------------------------------------
# ------------------ PyqtGraphWidget -------------------
# ------------------------------------------------------
# https://www.pythonguis.com/tutorials/pyside-plotting-pyqtgraph/

from PySide2.QtWidgets import QWidget, QVBoxLayout
from PySide2.QtGui import QFont
import pyqtgraph as pg  # pip install pyqtgraph


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

    def apply_styles(self):
        self.font_pen = pg.mkPen(color=self.font_color)

        self.title_font.setPixelSize(self.font_size_title)
        self.legend_font.setPixelSize(self.font_size_legend)
        self.axes_font.setPixelSize(self.font_size_axes)

        self.graphWidget.getAxis("bottom").setStyle(tickFont=self.axes_font)
        self.graphWidget.getAxis("left").setStyle(tickFont=self.axes_font)
        self.graphWidget.getAxis("bottom").setTextPen(self.font_pen)
        self.graphWidget.getAxis("left").setTextPen(self.font_pen)

        self.styles = {"color": "#000", "font-size": str(self.font_size_axes) + "px"}

        self.add_title()
        self.add_axes_titles()

        if self.legend_item:
            self.legend_item.setLabelTextSize(str(self.font_size_legend) + 'pt')

    def clear(self):
        self.graphWidget.clear()

    def add_legend(self):
        self.legend_item = self.graphWidget.addLegend()
        self.legend_item.setLabelTextSize(str(self.font_size_legend) + 'pt')
        self.legend_item.setLabelTextColor(pg.mkColor(0, 0, 0))

    def set_xticks(self, ticks) -> None:  # pragma: no cover
        self.graphWidget.getAxis("bottom").setTicks(ticks)

    def set_limits(self, x_min, x_max, y_min, y_max):  # pragma: no cover
        self.graphWidget.setXRange(x_min, x_max, padding=0)
        self.graphWidget.setYRange(y_min, y_max, padding=0)

    def enable_auto_range(self):  # pragma: no cover
        self.graphWidget.enableAutoRange()
