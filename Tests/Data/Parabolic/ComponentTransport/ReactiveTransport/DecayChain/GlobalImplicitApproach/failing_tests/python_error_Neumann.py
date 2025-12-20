# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class BC(OpenGeoSys.BoundaryCondition):
    def getFlux(self, _t, _coords, _primary_vars):
        msg = "this exception is thrown on purpose"
        raise RuntimeError(msg)


bc = BC()
