/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_LINALG_SPARSITYPATTERN_H
#define MATHLIB_LINALG_SPARSITYPATTERN_H

#include <vector>

#include "ProcessLib/NumericsConfig.h"

namespace MathLib
{
/// A vector telling how many nonzeros there are in each global matrix row.
using SparsityPattern = std::vector<GlobalIndexType>;
}

#endif // MATHLIB_LINALG_SPARSITYPATTERN_H
