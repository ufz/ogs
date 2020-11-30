# model of the evolving glacier extensions (length and height)
# parameterized, independent of concrete geometry

import numpy as np
import matplotlib.pyplot as plt

from math import pi, sin, cos, sinh, cosh, sqrt

gravity = 9.81  # m/s²


class glacier:
    def __init__(self, L_dom, L_max, H_max, x_0, t_0, t_1):
        self.rho_ice = 900  # kg/m³
        self.L_dom = L_dom
        self.L_max = L_max
        self.H_max = H_max
        self.x_0 = x_0
        self.t_0 = t_0
        self.t_1 = t_1

    def normalstress(self, x, t):
        return -self.rho_ice * gravity * self.local_height(x, t)

    # analytical function for the glacier's shape
    def local_height(self, x, t):
        l = self.length(t)
        if l == 0:
            return 0 * x
        else:
            xi = (x - self.x_0) / l
            xi = np.array(xi)
            xi[xi > 1] = 1.0
            return self.height(t) * np.sqrt(1 - xi ** 1)

    def height(self, t):
        return self.H_max * (t - self.t_0) / self.t_1

    def length(self, t):
        return self.L_max * (t - self.t_0) / self.t_1

    def printMaxLoads(self):
        print("Maximal normal stress due to glacier load: ")
        print(self.normalstress(0, self.t_1) / 1e6, "MPa")

    def plotEvolution(self):
        tRange = np.linspace(self.t_0, self.t_1, 11)
        xRange = np.linspace(self.x_0, self.x_0 + self.L_dom, 500)
        yRangeBefore = 0

        fig, ax = plt.subplots()
        ax.set_title("Glacier evolution")
        for t in tRange:
            yRange = self.local_height(xRange, t)
            ax.plot(xRange, yRange, label=t)
            ax.fill_between(xRange, yRangeBefore, yRange)
            yRangeBefore = yRange
        ax.set_xlabel("$x$ / m")
        ax.set_ylabel("height")
        ax.grid()
        fig.legend()
        fig.savefig("glacier.pdf")

        plt.show()

        fig, ax = plt.subplots()
        ax.plot(tRange, self.height(tRange))
        ax.set_xlabel("$t$ / a")
        ax.set_ylabel("height")
        ax.grid()
