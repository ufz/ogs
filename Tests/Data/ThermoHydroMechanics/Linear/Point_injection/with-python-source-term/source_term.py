# SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
# SPDX-License-Identifier: BSD-3-Clause

try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys


class SourceTerm(OpenGeoSys.SourceTerm):
    def getFlux(self, _t, coords, _primary_vars):
        _x, _y, _z = coords
        value = 150
        Jac = [0.0] * 4
        return (value, Jac)


source_term = SourceTerm()
