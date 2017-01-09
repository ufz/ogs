/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GlobalMatrixVectorTypes.h"

namespace MathLib
{

struct MatrixSpecifications
{
    MatrixSpecifications(std::size_t const nrows_, std::size_t const ncols_,
                         std::vector<GlobalIndexType> const*const ghost_indices_,
                         GlobalSparsityPattern const*const sparsity_pattern_)
        : nrows(nrows_), ncols(ncols_), ghost_indices(ghost_indices_)
        , sparsity_pattern(sparsity_pattern_)
    {
    }

    std::size_t const nrows;
    std::size_t const ncols;
    std::vector<GlobalIndexType> const*const ghost_indices;
    GlobalSparsityPattern const*const sparsity_pattern;
};

} // namespace MathLib
