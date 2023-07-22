import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Ellipse
from typing import List
import pandas as pd
from LineSegment import LineSegment
from ellipse import LsqEllipse


class Kepler:
    def __init__(self):
        self.r = 1
        self.Init()

        # Store intersecting points (for fitting ellispe later)
        self.ixs = []
        self.iys = []

        self.proc()
        self.fit()

    def Init(self):
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
        self.earthY = self.r * np.sin(theta)
        self.ax.plot(x, self.earthY, 'gray')

    def circle_intersection(self, line : LineSegment) -> List[float]:
        return [self.r * np.cos(line.angle), self.r * np.sin(line.angle)]


    def point(self, ang1, ang2) -> LineSegment:
        l1 = LineSegment(ang1, 0, 0, 2)
        l1.draw(self.ax, color='gray', alpha = 0.5)
        ix, iy = self.circle_intersection(l1)
        plt.plot(ix, iy, '.', color='gray')
        l2 = LineSegment(ang2, ix, iy, 2)
        l2.draw(self.ax, color='gray', alpha = 0.5)
        return l2

    def proc(self):
        df = pd.read_csv("data.csv")

        for index, row in df.iterrows():
            ang1s = np.array(df["HL Earth"])
            ang2s = np.array(df["GL Mars"])
        
        for i in range(len(ang1s)-1):
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
            edgecolor='b', fc='None', lw=2, label='Fit', zorder=2
        )
        self.ax.add_patch(ellipse)

        plt.xlabel('$X_1$')
        plt.ylabel('$X_2$')

        plt.legend()
        plt.show()

if __name__ == "__main__":
    Kepler()
