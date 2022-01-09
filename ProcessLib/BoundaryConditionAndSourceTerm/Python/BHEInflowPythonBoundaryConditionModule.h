/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <pybind11/pybind11.h>

#include "BaseLib/ExportSymbol.h"

namespace ProcessLib
{
//! Creates BHE Inflow Python bindings for the Python BC class.
OGS_EXPORT_SYMBOL void bheInflowpythonBindBoundaryCondition(
    pybind11::module& m);
}  // namespace ProcessLib
