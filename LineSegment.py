"""
    Class for represting a line segment
    Code By: V DHEERAJ SHENOY
"""

import matplotlib.pyplot as plt
from typing import List
import numpy as np

class LineSegment:
    """
        Takes the following argument:
            angle - in degrees
            x0 - starting x point
            y0 - starting y point
            length - length of the line segment
    """
    def __init__(self, angle, x0 = 0, y0 = 0, length = 1):
        angle %= 360
        self.angle = np.radians(angle)
        self.x0 = x0
        self.y0 = y0
        self.length = length
        self.x = self.x0 + self.length * np.cos(self.angle)
        self.y = self.y0 + self.length * np.sin(self.angle)
    
    def draw(self, ax, **kwargs):
        ax.plot([self.x0, self.x], [self.y0, self.y], **kwargs)

    def intersect(self, other : 'LineSegment'):
        x1,y1 = other.x0, other.y0
        x2,y2 = other.x, other.y
        x3,y3 = self.x0, self.y0
        x4,y4 = self.x, self.y

        denom = (y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)
        if denom == 0: # parallel
            return None
        ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3)) / denom
        if ua < 0 or ua > 1: # out of range
            return None
        ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3)) / denom
        if ub < 0 or ub > 1: # out of range
            return None
        x = x1 + ua * (x2-x1)
        y = y1 + ua * (y2-y1)
        return (x,y)
