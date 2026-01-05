# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

p_flux_in = 1e-2
p_0 = 1e5


# Source Terms
## Pressure
class BC_p_D(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, _t, _coords, _node_id, _primary_vars):
        return (True, p_0)


class BC_p_N(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, _coords, _primary_vars):
        Jac = [0.0, 0.0, 0.0]
        return (True, p_flux_in, Jac)


bc_p_D = BC_p_D()
bc_p_N = BC_p_N()
