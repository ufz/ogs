# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

from math import cos, cosh, pi, sin, sinh

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

a = 2.0 * pi / 3.0


# analytical solution used to set the Dirichlet BCs
def solution(x, y):
    return sin(a * x) * sinh(a * y)


# gradient of the analytical solution used to set the Neumann BCs
def grad_solution(x, y):
    return a * cos(a * x) * sinh(a * y), a * sin(a * x) * cosh(a * y)


# Dirichlet BCs
class BCTop(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        x, y, z = coords
        assert y == 1.0
        assert z == 0.0
        value = solution(x, y)
        return (True, value)


class BCLeft(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        x, y, z = coords
        assert x == 0.0
        assert z == 0.0
        value = solution(x, y)
        return (True, value)


class BCBottom(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, coords, _node_id, _primary_vars):
        x, y, z = coords
        assert y == 0.0
        assert z == 0.0
        value = solution(x, y)
        return (True, value)


# Neumann BC
class BCRight(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, coords, _primary_vars):
        x, y, z = coords
        assert x == 1.0
        assert z == 0.0
        value = grad_solution(x, y)[0]
        Jac = [0.0]  # value does not depend on primary variable
        return (True, value, Jac)


# instantiate BC objects referenced in OpenGeoSys' prj file
bc_top = BCTop()
bc_right = BCRight()
bc_bottom = BCBottom()
bc_left = BCLeft()
