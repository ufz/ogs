/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <pybind11/pybind11.h>

namespace ProcessLib
{
//! Creates Python bindings for the Python BC class.
void pythonBindBoundaryCondition(pybind11::module& m);
}  // namespace ProcessLib
