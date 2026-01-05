# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

rhoLR = 1000
cpL = 2000
v_x = 1.5e-6


class BC_Heat(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, _coords, primary_vars):
        temperature = primary_vars[2]
        print(temperature)
        value = -rhoLR * cpL * temperature * v_x
        return (True, value, [0.0, 0.0, 0.0, 0.0, 0.0])


# instantiate the BC objects used by OpenGeoSys
neumann_correction = BC_Heat()
