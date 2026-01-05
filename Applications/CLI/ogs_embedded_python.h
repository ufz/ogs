// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <pybind11/embed.h>

#include "BaseLib/ExportSymbol.h"

namespace ApplicationsLib
{
/// Sets up an embedded Python interpreter and makes sure that the OpenGeoSys
/// Python module is not removed by the linker.
OGS_EXPORT_SYMBOL pybind11::scoped_interpreter setupEmbeddedPython();

/// Checks for activated and matching virtual environment and adds it to
/// sys.path.
OGS_EXPORT_SYMBOL void setupEmbeddedPythonVenvPaths();

}  // namespace ApplicationsLib
