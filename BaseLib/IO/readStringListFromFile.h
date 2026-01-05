// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string>
#include <vector>

namespace BaseLib
{
namespace IO
{
/// Reads non-empty lines from a list of strings from a file into a vector
std::vector<std::string> readStringListFromFile(std::string const& filename);
}
}  // namespace BaseLib
