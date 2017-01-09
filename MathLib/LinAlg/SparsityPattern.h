/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_LINALG_SPARSITYPATTERN_H
#define MATHLIB_LINALG_SPARSITYPATTERN_H

#include <vector>

namespace MathLib
{
/// A vector telling how many nonzeros there are in each global matrix row.
template <typename IndexType>
using SparsityPattern = std::vector<IndexType>;
}

#endif // MATHLIB_LINALG_SPARSITYPATTERN_H
