# model of the ideal gas inside the chamber

from math import pi

R_piston = 0.05

h_chamber_0 = 0.1
V_chamber_0 = pi * R_piston**2 * h_chamber_0

p_chamber_0 = 1e5


def p_chamber(displacement_y):
    V_chamber = pi * R_piston**2 * (h_chamber_0 + displacement_y)
    return p_chamber_0 * V_chamber_0 / V_chamber

def dp_chamber_du_y(displacement_y):
    return -p_chamber_0 * V_chamber_0 / pi / R_piston**2 / \
            (h_chamber_0 + displacement_y)**2
