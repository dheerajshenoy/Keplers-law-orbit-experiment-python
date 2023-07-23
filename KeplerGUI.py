#!/usr/bin/env python3
"""
    The program is based on a lab experiment from the link in the README.md. The procedure is simple:
    There are two data given, Heliocentric Longitude of Earth (HL) and Geocentric Longitude of Mars (GL).
    Both these data are in degrees. We first draw a circle of a radius depicting the Earth's orbit.
    Then we mark the center of this circle as the sun. Next, we find the intersection of the line that
    makes an angle (HL) from the sun and move to this intersecting point. We draw a line making an angle
    GL from this point. We go to the next data and repeat the same procedure. We find that these two lines
    having angles (GL1, GL2) intersect at a point. We repeat this procedure for the remaining data.
    Next, we find the best fit of an ellipse through these set of intersecting points. This gives the
    orbit of Mars(in this case).

    Code By: V DHEERAJ SHENOY
"""

from PyQt6.QtGui import QColor
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse
from typing import List
import pandas as pd
from LineSegment import LineSegment
from ellipse import LsqEllipse
import sys
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QGridLayout,
                             QGroupBox, QMenuBar, QPushButton, QLineEdit,
                             QColorDialog, QLabel, QSizePolicy, QSlider,
                             QSplitter, QMainWindow, QApplication, QWidget,
                             QVBoxLayout, QHBoxLayout)
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class Kepler:
    def __init__(self, fig, ax, filename = None):
        self.fig = fig
        self.ax = ax
        self.color_sun = 'y'
        self.color_earth_orbit = 'k'
        self.shape_line_circle_intersection = 'o'
        self.shape_intersection_point = 'o'
        self.color_line_circle_intersection = 'gray'
        self.color_line_from_sun = 'gray'
        self.color_line_from_earth = 'green'
        self.color_planet_orbit = 'gold'
        self.color_intersection_point = 'purple'
        self.alpha_earth_orbit = 0.2
        self.alpha_line_from_sun = 0.2
        self.alpha_line_from_earth = 0.2
        self.alpha_intersection = 0.2
        self.alpha_line_circle_intersection = 0.2
        self.alpha_planet_orbit = 0.2
        self.r = 1 # radius of earth's orbit
        self.Init()
        self.l1s = []
        self.l2s = []
        self.lineCircleIntersectionList = []
        self.intersectionList = []

        # Store intersecting points (for fitting ellispe later)
        self.ixs = []
        self.iys = []

        if filename:
            self.proc(filename) # read data from the csv file and mark the intersecting points
            self.fit() # Best fit an ellipse through these intersecting points

    def DataFile(self, file):
        self.proc(file)
        self.fit()

    def Init(self):
        """
            Initialising function that draws the sun at the center, draws Earth's orbit
        """
        self.ax.set_xlim(-2, 2)
        self.ax.set_ylim(-2, 2)
        self.ax.set_aspect('equal')

        # Draw Sun
        self.sun = Circle(xy=(0, 0), radius = 0.1, color=self.color_sun)
        self.ax.add_artist(self.sun)

        # Draw Earth's orbit
        self.earth_orbit()

    def earth_orbit(self):
        theta = np.linspace(0, 2 * np.pi, 1000)
        x = self.r * np.cos(theta)
        y = self.r * np.sin(theta)
        self.earthorbit = self.ax.plot(x, y, color=self.color_earth_orbit,
                                       label = 'Earth\'s Orbit', alpha = self.alpha_earth_orbit)
        type(self.earthorbit)

    def circle_intersection(self, line : LineSegment) -> List[float]:
        # Function that returns the coordinates (x, y) of the point of intersection of the HL line with the Earth's orbit
        return [self.r * np.cos(line.angle), self.r * np.sin(line.angle)]

    def point(self, ang1, ang2) -> LineSegment:
        """
            l1 is the line drawn from the sun { center (0, 0) }.
            l2 is the line drawn from the point of intersection of the l1 line with the Earth's orbit
            Function returns l2 line
        """
        l1 = LineSegment(ang1, 0, 0, 2)
        self.l1s.append(l1.draw(self.ax, color=self.color_line_from_sun, alpha = self.alpha_line_from_sun))

        ix, iy = self.circle_intersection(l1)
        self.lineCircleIntersectionList.append(self.ax.plot(ix, iy,
                                                            marker=self.shape_line_circle_intersection,
                                                            color=self.color_line_circle_intersection,
                                                            alpha = self.alpha_line_circle_intersection))
        l2 = LineSegment(ang2, ix, iy, 2)
        self.l2s.append(l2.draw(self.ax, color=self.color_line_from_earth, alpha = self.alpha_line_from_earth))
        return l2

    def proc(self, csv_filename):
        """
            Reads the csv data from `csv_filename` file and reads the angles and runs the procedure
            of plotting the lines and intersecting points of these lines and stores all these
            intersecting points.
        """
        df = pd.read_csv(csv_filename)

        for index, row in df.iterrows():
            ang1s = np.array(df["HL Earth"])
            ang2s = np.array(df["GL Mars"])

        for i in range(len(ang1s)-1):
            # We need two set of angles to find the intersection of the lines
            ang1 = ang1s[i]
            ang2 = ang2s[i]
            ang3 = ang1s[i + 1]
            ang4 = ang2s[i + 1]
            l1 = self.point(ang1, ang2)
            l2 = self.point(ang3, ang4)
            ipt = l1.intersect(l2)
            if ipt is not None:
                self.intersectionList.append(self.ax.plot(*ipt, self.color_intersection_point,
                                                          marker=self.shape_intersection_point, markersize=10))
                self.ixs.append(ipt[0])
                self.iys.append(ipt[1])

    def fit(self):
        x = np.array(self.ixs)
        y = np.array(self.iys)
        X = np.array(list(zip(x, y)))
        reg = LsqEllipse().fit(X)
        center, width, height, phi = reg.as_parameters()

        # print(f'center: {center[0]:.3f}, {center[1]:.3f}')
        # print(f'width: {width:.3f}')
        # print(f'height: {height:.3f}')
        # print(f'phi: {phi:.3f}')
        #
        self.planetOrbit = Ellipse(
                xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
                edgecolor=self.color_planet_orbit, fc='None', lw=2, label='Mar\'s Orbit', zorder=2,
                alpha = self.alpha_planet_orbit)
        self.ax.add_patch(self.planetOrbit)

    def change_earth_orbit_color(self, color):
        self.color_earth_orbit = color
        self.earthorbit[0].set_color(color)
        self.fig.canvas.draw()

    def change_intersection_point_color(self, color):
        self.color_intersection_point = color
        for i in self.intersectionList:
            i[0].set_color(color)
        self.fig.canvas.draw()

    def change_line_circle_intersection_color(self, color):
        self.color_line_circle_intersection = color
        for i in self.lineCircleIntersectionList:
            i[0].set_color(color)
        self.fig.canvas.draw()

    def change_line_from_sun_color(self, color):
        self.color_sun = color
        for i in self.l1s:
            i[0].set_color(color)
        self.fig.canvas.draw()

    def change_planet_orbit_color(self, color):
        self.color_planet_orbit = color

        try:
            self.planetOrbit.set_edgecolor(color)
        except AttributeError:
            pass
        self.fig.canvas.draw()

    def change_line_from_earth_color(self, color):
        self.color_line_from_earth = color
        for i in self.l2s:
            i[0].set_color(color)
        self.fig.canvas.draw()

    def change_intersection_point_shape(self, shape):
        self.shape_intersection_point = shape
        for i in self.intersectionList:
            i[0].set_marker(shape)
        self.fig.canvas.draw()

    def change_circle_intersection_point_shape(self, shape):
        self.shape_line_circle_intersection = shape
        for i in self.lineCircleIntersectionList:
            i[0].set_marker(shape)
        self.fig.canvas.draw()

    def change_earth_orbit_alpha(self, value):
        value = value / 100
        self.alpha_earth_orbit = value
        self.earthorbit[0].set_alpha(value)
        self.fig.canvas.draw()

    def change_planet_orbit_alpha(self, value):
        value = value / 100
        self.alpha_planet_orbit = value
        try:
            self.planetOrbit.set_alpha(value)
            self.fig.canvas.draw()
        except AttributeError:
            pass

    def change_earth_to_planet_line_alpha(self, value):
        value = value / 100
        self.alpha_line_from_earth = value
        for i in self.l2s:
            i[0].set_alpha(value)
        self.fig.canvas.draw()


    def change_sun_to_earth_line_alpha(self, value):
        value = value / 100
        self.alpha_line_from_sun = value
        for i in self.l1s:
            i[0].set_alpha(value)
        self.fig.canvas.draw()

    def change_intersection_point_alpha(self, value):
        value = value / 100
        self.alpha_intersection = value
        for i in self.intersectionList:
            i[0].set_alpha(value)
        self.fig.canvas.draw()

    def change_line_circle_intersection_alpha(self, value):
        value = value / 100
        self.alpha_line_circle_intersection = value
        for i in self.lineCircleIntersectionList:
            i[0].set_alpha(value)
        self.fig.canvas.draw()

class ColorBar(QWidget):
    def __init__(self, text, default_color, **kwargs):
        super().__init__(**kwargs)
        self.layout = QHBoxLayout()
        self.setLayout(self.layout)
        self.btn = QPushButton(text)
        self.btn.clicked.connect(self.getColor)
        self.color = default_color
        self.btn.setStyleSheet("background-color: {}".format(self.color))
        self.layout.addWidget(self.btn, stretch=True)
        self.dialog = QColorDialog()

    def getColor(self):
        self.color = self.dialog.getColor().name()
        if self.function is not None:
            self.function(self.color)

        color = QColor(self.color)
        print(color.lightness())
        if color.lightness() >= 20 and color.lightness() <= 255:
            self.btn.setStyleSheet("color: {}; background-color: {}".format(color.lighter().name(), self.color))
        elif color.lightness() < 20:
            self.btn.setStyleSheet("color: {}; background-color: {}".format(color.darker().name(), self.color))

    def color(self):
        return self.color

    def connect(self, function):
        if function is not None:
            self.function = function

class ComboBar(QWidget):
    def __init__(self, text, items : List[str] = [], default = None, **kwargs):
        super().__init__(**kwargs)

        self.layout = QHBoxLayout()
        self.setLayout(self.layout)

        self.label = QLabel(text)
        self.box = QComboBox()

        if len(items) > 0:
            self.box.addItems(items)

        if default:
            self.box.setCurrentIndex(items.index(default))
        self.layout.addWidget(self.label)
        self.layout.addWidget(self.box)

    def connect(self, function):
        if function is not None:
            self.box.currentTextChanged.connect(function)

class AlphaBar(QWidget):
    def __init__(self, text, From, To, Interval, **kwargs):
        super().__init__(**kwargs)

        self.layout = QHBoxLayout()

        self.label = QLabel(text)

        self.slider =QSlider(Qt.Orientation.Horizontal, )
        self.slider.setSingleStep(25)
        self.slider.setTickPosition(self.slider.TickPosition.TicksAbove)
        self.slider.setMinimum(From)
        self.slider.setMaximum(To)
        self.slider.setTickInterval(Interval)

        self.setLayout(self.layout)

        self.layout.addWidget(self.label)
        self.layout.addWidget(self.slider)

    def connect(self, function):
        if function is not None:
            self.slider.valueChanged.connect(function)

class SideBar(QWidget):
    def __init__(self, win : QMainWindow, kepler : Kepler, **kwargs):
        super().__init__(**kwargs)

        self.win = win
        self.kepler = kepler
        self.setMinimumWidth(400)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # [Data File Groupbox]
        self.gbx_data = QGroupBox("Data")
        self.layout.addWidget(self.gbx_data)
        self.gbx_data_layout = QGridLayout()
        self.gbx_data.setLayout(self.gbx_data_layout)

        self.file_label = QLabel("File")
        self.file_browse_button = QPushButton("Browse...")
        self.file_browse_button.clicked.connect(self.data_file_func)

        self.fileWidget = QWidget()
        self.fileLayout = QHBoxLayout()
        self.fileLayout.addWidget(self.file_label)
        self.fileLayout.addWidget(self.file_browse_button, stretch=True)

        self.fileWidget.setLayout(self.fileLayout)
        self.gbx_data_layout.addWidget(self.fileWidget)

        # [Toggle Groupbox]
        self.gbx_toggle = QGroupBox("Toggle")
        self.layout.addWidget(self.gbx_toggle)
        self.gbx_toggle_layout = QGridLayout()
        self.gbx_toggle.setLayout(self.gbx_toggle_layout)

        self.cbx_grid = QCheckBox("Grid")
        self.cbx_earth_orbit = QCheckBox("Earth Orbit")
        self.cbx_sun = QCheckBox("Sun")
        self.cbx_planet_orbit = QCheckBox("Planet Orbit")
        self.cbx_earth_orbit_intersection_point = QCheckBox("Earth Orbit Intersection Points")
        self.cbx_line_from_sun = QCheckBox("Line from sun")
        self.cbx_lines = QCheckBox("Lines from Earth Orbit")
        self.cbx_intersection_point = QCheckBox("Line Intersection Points")
        self.cbx_legend = QCheckBox("Legend")
        self.cbx_axis = QCheckBox("Axis")
        self.cbx_toolbar = QCheckBox("Toolbar")

        self.gbx_toggle_layout.addWidget(self.cbx_grid, 0, 0)
        self.gbx_toggle_layout.addWidget(self.cbx_earth_orbit, 0, 1)
        self.gbx_toggle_layout.addWidget(self.cbx_sun, 1, 0)
        self.gbx_toggle_layout.addWidget(self.cbx_planet_orbit, 1, 1)
        self.gbx_toggle_layout.addWidget(self.cbx_earth_orbit_intersection_point, 2, 0)
        self.gbx_toggle_layout.addWidget(self.cbx_line_from_sun, 2, 1)
        self.gbx_toggle_layout.addWidget(self.cbx_lines, 3, 0)
        self.gbx_toggle_layout.addWidget(self.cbx_toolbar, 3, 1)
        self.gbx_toggle_layout.addWidget(self.cbx_intersection_point, 4, 0)
        self.gbx_toggle_layout.addWidget(self.cbx_legend, 4, 1)
        self.gbx_toggle_layout.addWidget(self.cbx_axis, 5, 0)

        # [Shape Groupbox]
        self.gbx_shape = QGroupBox("Shape")
        self.layout.addWidget(self.gbx_shape)

        self.gbx_shape_layout = QVBoxLayout()
        self.gbx_shape.setLayout(self.gbx_shape_layout)

        # Get list of all markers
        # self.Markers = []
        # for k in MarkerStyle.markers.keys():
        #     self.Markers.append(str(k))

        self.Markers = [".", "o", "x", "X", ",", "v", "^", "<", ">", "1",
                        "2", "3", "4", "8", "s", "p", "P", "*", "h", "H",
                        "+", "D", "d", "|", "_"]

        self.circle_intersection_point_shape = ComboBar("Circle Intersection", self.Markers, '.')
        self.circle_intersection_point_shape.connect(self.kepler.change_circle_intersection_point_shape)

        self.intersection_point_shape = ComboBar("Intersection", self.Markers, '.')
        self.intersection_point_shape.connect(self.kepler.change_intersection_point_shape)

        self.gbx_shape_layout.addWidget(self.circle_intersection_point_shape)
        self.gbx_shape_layout.addWidget(self.intersection_point_shape)

        # [Color Groupbox]
        self.gbx_color = QGroupBox("Color")
        self.layout.addWidget(self.gbx_color)

        self.gbx_color_layout = QGridLayout()
        self.gbx_color.setLayout(self.gbx_color_layout)

        self.earth_orbit_color = ColorBar("Earth Orbit", self.kepler.color_earth_orbit)
        self.earth_orbit_color.connect(self.kepler.change_earth_orbit_color)

        self.planet_orbit_color = ColorBar("Planet Orbit", self.kepler.color_planet_orbit)
        self.planet_orbit_color.connect(self.kepler.change_planet_orbit_color)

        self.circle_intersection_point_color = ColorBar("Point", self.kepler.color_line_circle_intersection)
        self.circle_intersection_point_color.connect(self.kepler.change_line_circle_intersection_color)

        self.intersection_point_color = ColorBar("Intersection", self.kepler.color_intersection_point)
        self.intersection_point_color.connect(self.kepler.change_intersection_point_color)

        self.line_from_sun_color = ColorBar("Line from sun", self.kepler.color_line_from_sun)
        self.line_from_sun_color.connect(self.kepler.change_line_from_sun_color)

        self.line_to_planet_color = ColorBar("Line To Planet", self.kepler.color_line_from_earth)
        self.line_to_planet_color.connect(self.kepler.change_line_from_earth_color)

        self.gbx_color_layout.addWidget(self.earth_orbit_color, 0, 0)
        self.gbx_color_layout.addWidget(self.planet_orbit_color, 0, 1)
        self.gbx_color_layout.addWidget(self.circle_intersection_point_color, 1, 0)
        self.gbx_color_layout.addWidget(self.intersection_point_color, 1, 1)
        self.gbx_color_layout.addWidget(self.line_from_sun_color, 2, 0)
        self.gbx_color_layout.addWidget(self.line_to_planet_color, 2, 1)

        self.setupBools()

        # [Alpha Groupbox]

        self.gbx_alpha = QGroupBox("alpha")

        self.layout.addWidget(self.gbx_alpha)
        self.gbx_alpha_layout = QGridLayout()
        self.gbx_alpha.setLayout(self.gbx_alpha_layout)

        self.sl_earth_orbit_alpha = AlphaBar("Earth Orbit", 0, 100, 25)
        self.sl_earth_orbit_alpha.connect(self.kepler.change_earth_orbit_alpha)

        self.sl_planet_orbit_alpha = AlphaBar("Planet Orbit", 0, 100, 25)
        self.sl_planet_orbit_alpha.connect(self.kepler.change_planet_orbit_alpha)

        self.sl_earth_to_planet_line_alpha = AlphaBar("Line from Earth", 0, 100, 25)
        self.sl_earth_to_planet_line_alpha.connect(self.kepler.change_earth_to_planet_line_alpha)

        self.sl_sun_to_earth_line_alpha = AlphaBar("Line from Sun", 0, 100, 25)
        self.sl_sun_to_earth_line_alpha.connect(self.kepler.change_sun_to_earth_line_alpha)

        self.sl_intersection_alpha = AlphaBar("Intersection", 0, 100, 25)
        self.sl_intersection_alpha.connect(self.kepler.change_intersection_point_alpha)

        self.sl_circle_line_intersection_alpha = AlphaBar("Circle Intersection", 0, 100, 25)
        self.sl_circle_line_intersection_alpha.connect(self.kepler.change_line_circle_intersection_alpha)

        self.gbx_alpha_layout.addWidget(self.sl_earth_orbit_alpha, 0, 0)
        self.gbx_alpha_layout.addWidget(self.sl_planet_orbit_alpha, 0, 1)
        self.gbx_alpha_layout.addWidget(self.sl_intersection_alpha, 1, 0)
        self.gbx_alpha_layout.addWidget(self.sl_sun_to_earth_line_alpha, 1, 1)
        self.gbx_alpha_layout.addWidget(self.sl_earth_to_planet_line_alpha, 2, 0)
        self.gbx_alpha_layout.addWidget(self.sl_circle_line_intersection_alpha, 2, 1)

    def updatecanvas(self):
        self.kepler.fig.canvas.draw()

    def setupBools(self):
        self.b_grid = True
        self.b_toolbar = True
        self.b_sun = True
        self.b_earth_orbit = True
        self.b_planet_orbit = True
        self.b_circle_intersection_points = True
        self.b_intersection_points = True
        self.b_earth_orbit_intersection = True
        self.b_line_from_sun = True
        self.b_line_from_earth_orbit = True
        self.b_legend = True
        self.b_axis = False

        self.cbx_grid.clicked.connect(self.toggle_grid)
        self.cbx_earth_orbit.clicked.connect(self.toggle_earth_orbit)
        self.cbx_sun.clicked.connect(self.toggle_sun)
        self.cbx_planet_orbit.clicked.connect(self.toggle_planet_orbit)
        self.cbx_earth_orbit_intersection_point.clicked.connect(self.toggle_earth_orbit_intersection_point)
        self.cbx_line_from_sun.clicked.connect(self.toggle_line_from_sun)
        self.cbx_lines.clicked.connect(self.toggle_line_from_earth)
        self.cbx_intersection_point.clicked.connect(self.toggle_intersection_points)
        self.cbx_legend.clicked.connect(self.toggle_legend)
        self.cbx_toolbar.clicked.connect(self.toggle_toolbar)
        self.cbx_axis.clicked.connect(self.toggle_axis)

        self.cbx_grid.setChecked(self.b_grid)
        self.toggle_grid(self.b_grid)

        self.cbx_earth_orbit.setChecked(self.b_earth_orbit)

        self.cbx_sun.setChecked(self.b_sun)

        self.cbx_planet_orbit.setChecked(self.b_planet_orbit)
        self.cbx_planet_orbit.setEnabled(False)

        self.cbx_earth_orbit_intersection_point.setChecked(self.b_earth_orbit_intersection)
        self.cbx_earth_orbit_intersection_point.setEnabled(False)

        self.cbx_line_from_sun.setChecked(self.b_line_from_sun)
        self.cbx_line_from_sun.setEnabled(False)

        self.cbx_lines.setChecked(self.b_line_from_earth_orbit)

        self.cbx_toolbar.setChecked(self.b_toolbar)

        self.cbx_intersection_point.setChecked(self.b_intersection_points)
        self.cbx_intersection_point.setEnabled(False)

        self.cbx_legend.setChecked(self.b_legend)
        self.toggle_legend(self.b_legend)

        self.cbx_axis.setChecked(self.b_axis)
        self.toggle_axis(self.b_axis)

    def data_file_func(self):
        filepath = QFileDialog().getOpenFileName(self, "Open File", "~", "CSV Files (*.csv)")
        data_file = filepath[0]
        self.file_browse_button.setText(data_file.split("/")[-1])
        self.kepler.DataFile(data_file)
        self.cbx_planet_orbit.setEnabled(True)
        self.cbx_line_from_sun.setEnabled(True)
        self.cbx_intersection_point.setEnabled(True)
        self.cbx_earth_orbit_intersection_point.setEnabled(True)
        self.updatecanvas()

    def toggle_grid(self, value):
        self.b_grid = value
        self.kepler.ax.grid(value)
        self.updatecanvas()

    def toggle_earth_orbit(self, value):
        self.b_earth_orbit = value

        if value:
            self.kepler.earth_orbit()
        else:
            j = self.kepler.earthorbit.pop(0)
            j.remove()
        self.updatecanvas()

    def toggle_sun(self, value):
        self.b_sun = value
        self.kepler.sun.set_visible(value)
        self.updatecanvas()

    def toggle_planet_orbit(self, value):
        self.b_planet_orbit = value
        self.kepler.planetOrbit.set_visible(value)
        self.updatecanvas()

    def toggle_earth_orbit_intersection_point(self, value):
        self.b_earth_orbit_intersection = value
        if len(self.kepler.lineCircleIntersectionList) > 0:
            if not value:
                for ic in self.kepler.lineCircleIntersectionList:
                    ic[0].set_visible(False)
            else:
                for ic in self.kepler.lineCircleIntersectionList:
                    ic[0].set_visible(True)
            self.updatecanvas()

    def toggle_line_from_sun(self, value):
        self.b_line_from_sun = value

        if len(self.kepler.l1s) > 0:
            if not value:
                for l1 in self.kepler.l1s:
                    l1[0].set_visible(False)
            else:
                for l1 in self.kepler.l1s:
                    l1[0].set_visible(True)

            self.updatecanvas()

    def toggle_line_from_earth(self, value):
        self.b_line_from_earth_orbit = value
        if len(self.kepler.l1s) > 0:
            if not value:
                for l2 in self.kepler.l2s:
                    l2[0].set_visible(False)
            else:
                for l2 in self.kepler.l2s:
                    l2[0].set_visible(True)

            self.updatecanvas()

    def toggle_intersection_points(self, value):
        self.b_intersection_points = value
        if len(self.kepler.intersectionList) > 0:
            if not value:
                for i in self.kepler.intersectionList:
                    i[0].set_visible(False)
            else:
                for i in self.kepler.intersectionList:
                    i[0].set_visible(True)

            self.updatecanvas()

    def toggle_legend(self, value):
        self.b_legend = value
        if value:
            self.kepler.ax.legend()
        else:
            self.kepler.ax.get_legend().remove()
        self.updatecanvas()

    def toggle_axis(self, value):
        self.b_axis_tick_labels = value
        self.kepler.ax.get_xaxis().set_visible(value)
        self.kepler.ax.get_yaxis().set_visible(value)

        if not value:
            self.cbx_grid.setEnabled(value)
        else:
            self.cbx_grid.setEnabled(True)
        self.updatecanvas()

    def toggle_toolbar(self, value):
        self.b_toolbar = value
        self.win.toolbar.setVisible(value)

class MainApplication(QMainWindow):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.fig, self.ax = plt.subplots()
        self.kepler = Kepler(self.fig, self.ax)
        self.InitGUI(self.kepler)
        self.show()

    def InitGUI(self, kepler : Kepler):
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.mainWidget = QWidget()
        self.mainLayout = QVBoxLayout()

        self.InitMenubar()
        
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.setCentralWidget(self.mainWidget)
        self.mainWidget.setLayout(self.mainLayout)

        self.sideBar = SideBar(self, kepler)
        self.canvasWidget = QWidget()
        self.canvasWidgetLayout = QVBoxLayout()
        self.canvasWidget.setLayout(self.canvasWidgetLayout)
        self.canvasWidgetLayout.addWidget(self.toolbar)
        self.canvasWidgetLayout.addWidget(self.canvas)

        self.splitSection = QSplitter(Qt.Orientation.Horizontal)

        self.splitSection.setStyleSheet('QSplitter::handle { width: 1px; background: gray; }')
        self.splitSection.setContentsMargins(0, 0, 0, 0)
        self.splitSection.addWidget(self.sideBar)
        self.splitSection.addWidget(self.canvasWidget)
        self.mainLayout.setContentsMargins(0, 0, 0, 0)
        self.mainLayout.addWidget(self.menubar)
        self.mainLayout.addWidget(self.splitSection)

    def InitMenubar(self):
        self.menubar = QMenuBar()
        self.menubar.addMenu("File")
        self.menubar.addMenu("Help")
        self.menubar.addMenu("Exit")

if __name__ == "__main__":
    app = QApplication(sys.argv)

    win = MainApplication()

    sys.exit(app.exec())
