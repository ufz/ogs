/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
