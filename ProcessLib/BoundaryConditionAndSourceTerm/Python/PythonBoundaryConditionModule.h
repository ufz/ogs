// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <pybind11/pybind11.h>

#include "BaseLib/ExportSymbol.h"

namespace ProcessLib
{
//! Creates Python bindings for the Python BC class.
OGS_EXPORT_SYMBOL void pythonBindBoundaryCondition(pybind11::module& m);
}  // namespace ProcessLib
