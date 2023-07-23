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

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse
from typing import List
import pandas as pd
from LineSegment import LineSegment
from ellipse import LsqEllipse

class Kepler:
    def __init__(self):
        self.r = 1 # radius of earth's orbit
        self.Init()

        # Store intersecting points (for fitting ellispe later)
        self.ixs = []
        self.iys = []

        self.proc("data.csv") # read data from the csv file and mark the intersecting points
        self.fit() # Best fit an ellipse through these intersecting points


    def Init(self):
        """
            Initialising function that draws the sun at the center, draws Earth's orbit
        """
        self.fig, self.ax = plt.subplots()
        self.ax.set_xlim(-2, 2)
        self.ax.set_ylim(-2, 2)
        self.ax.set_aspect('equal')
        self.ax.grid()

        # Draw Sun
        sun = Circle(xy=(0, 0), radius = 0.1, color='y')
        self.ax.add_artist(sun)

        # Draw Earth's orbit
        theta = np.linspace(0, 2 * np.pi, 1000)
        x = self.r * np.cos(theta)
        y = self.r * np.sin(theta)
        self.ax.plot(x, y, color='y', label = 'Earth\'s Orbit')

    def circle_intersection(self, line : LineSegment) -> List[float]:
        """
            Function that returns the coordinates (x, y) of the point of intersection of the HL line with the Earth's orbit
        """
        return [self.r * np.cos(line.angle), self.r * np.sin(line.angle)]

    def point(self, ang1, ang2) -> LineSegment:
        """
            l1 is the line drawn from the sun { center (0, 0) }.
            l2 is the line drawn from the point of intersection of the l1 line with the Earth's orbit
            Function returns l2 line
        """
        l1 = LineSegment(ang1, 0, 0, 2)
        l1.draw(self.ax, color='gray', alpha = 0.2)
        ix, iy = self.circle_intersection(l1)
        plt.plot(ix, iy, '.', color='gray')
        l2 = LineSegment(ang2, ix, iy, 2)
        l2.draw(self.ax, color='gray', alpha = 0.2)
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
                self.ax.plot(*ipt, 'k', marker='o', markersize=10)
                self.ixs.append(ipt[0])
                self.iys.append(ipt[1])

    def fit(self):
        x = np.array(self.ixs)
        y = np.array(self.iys)
        X = np.array(list(zip(x, y)))
        reg = LsqEllipse().fit(X)
        center, width, height, phi = reg.as_parameters()

        print(f'center: {center[0]:.3f}, {center[1]:.3f}')
        print(f'width: {width:.3f}')
        print(f'height: {height:.3f}')
        print(f'phi: {phi:.3f}')

        ellipse = Ellipse(
            xy=center, width=2*width, height=2*height, angle=np.rad2deg(phi),
            edgecolor='b', fc='None', lw=2, label='Mar\'s Orbit', zorder=2
        )
        self.ax.add_patch(ellipse)

        plt.legend()
        plt.show()

if __name__ == "__main__":
    Kepler()
