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

}  // namespace ApplicationsLib
