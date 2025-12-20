# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys
from math import atan2, cos, pi, sin, sqrt

a = 2.0  # Length
b = 1.0  # Height
x_tip = a / 4.0  # initial crack location
y_tip = b / 2.0  # initial crack location


E = 210.0e3  # MPa
Gc = 2.7
nu = 0.3
v = 1.5  # mm/s

kappa = (3.0 - nu) / (1.0 + nu)
mu = E / (2.0 * (1.0 + nu))

PlaneStress = False
if not PlaneStress:  # plane strain
    KI = float(sqrt(E * Gc / (1.0 - nu**2)))
else:  # plane stress
    KI = float(sqrt(E * Gc))


# Asymptotic solution used to set the Dirichlet BCS, X-dir
def AsymptoticSolution_X(x, y, t):
    return (
        KI
        / (2 * mu)
        * sqrt(
            sqrt((x - x_tip - v * t) * (x - x_tip - v * t) + (y - y_tip) * (y - y_tip))
            / (2 * pi)
        )
        * (kappa - cos(atan2(y - y_tip, x - x_tip - v * t)))
        * cos(atan2(y - y_tip, x - x_tip - v * t) / 2)
    )


# Asymptotic solution used to set the Dirichlet BCS, Y-dir
def AsymptoticSolution_Y(x, y, t):
    return (
        KI
        / (2 * mu)
        * sqrt(
            sqrt((x - x_tip - v * t) * (x - x_tip - v * t) + (y - y_tip) * (y - y_tip))
            / (2 * pi)
        )
        * (kappa - cos(atan2(y - y_tip, x - x_tip - v * t)))
        * sin(atan2(y - y_tip, x - x_tip - v * t) / 2)
    )


# Dirichlet BCs X-dir
class BCTop_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert y == b
        value = AsymptoticSolution_X(x, y, t)
        return (True, value)


class BCLeft_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert x == 0.0
        value = AsymptoticSolution_X(x, y, t)
        return (True, value)


class BCBottom_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert y == 0.0
        value = AsymptoticSolution_X(x, y, t)
        return (True, value)


class BCRight_X(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert x == a
        value = AsymptoticSolution_X(x, y, t)
        return (True, value)


# Dirichlet BCs Y-dir
class BCTop_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert y == b
        value = AsymptoticSolution_Y(x, y, t)
        return (True, value)


class BCLeft_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert x == 0.0
        value = AsymptoticSolution_Y(x, y, t)
        return (True, value)


class BCBottom_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert y == 0.0
        value = AsymptoticSolution_Y(x, y, t)
        return (True, value)


class BCRight_Y(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, _node_id, _primary_vars):
        x, y, _z = coords
        assert x == a
        value = AsymptoticSolution_Y(x, y, t)
        return (True, value)


# instantiate BC objects referenced in OpenGeoSys' prj file
bc_top_X = BCTop_X()
bc_right_X = BCRight_X()
bc_bottom_X = BCBottom_X()
bc_left_X = BCLeft_X()


bc_top_Y = BCTop_Y()
bc_right_Y = BCRight_Y()
bc_bottom_Y = BCBottom_Y()
bc_left_Y = BCLeft_Y()
