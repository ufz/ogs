// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
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
