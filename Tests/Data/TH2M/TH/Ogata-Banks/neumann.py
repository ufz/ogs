# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

rhoLR = 1000
cpL = 2000
v_x = 1.5e-6


# Component indices are positional: primary_vars[2] reads temperature in the
# TH2M monolithic component order (gas pressure, capillary pressure,
# temperature, displacements), and the Jacobian entry must sit at that same
# index of the returned dFlux list.
class BC_Heat(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, _coords, primary_vars):
        temperature = primary_vars[2]
        value = -rhoLR * cpL * temperature * v_x
        dvalue_dT = -rhoLR * cpL * v_x
        return (True, value, [0.0, 0.0, dvalue_dT, 0.0, 0.0])


# instantiate the BC objects used by OpenGeoSys
neumann_correction = BC_Heat()
