/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace MathLib
{
/// A vector telling how many nonzeros there are in each global matrix row.
template <typename IndexType>
using SparsityPattern = std::vector<IndexType>;
}
