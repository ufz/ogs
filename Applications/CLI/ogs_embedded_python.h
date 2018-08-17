/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <pybind11/embed.h>

namespace ApplicationsLib
{
/// Sets up an embedded Python interpreter and makes sure that the OpenGeoSys
/// Python module is not removed by the linker.
pybind11::scoped_interpreter setupEmbeddedPython();

}  // namespace ApplicationsLib
