/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
namespace SourceTerms
{
namespace Python
{
//! Creates Python bindings for the Python source term class.
OGS_EXPORT_SYMBOL void pythonBindSourceTerm(pybind11::module& m);
}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
