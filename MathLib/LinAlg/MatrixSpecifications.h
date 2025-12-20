// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "GlobalMatrixVectorTypes.h"

namespace MathLib
{

struct MatrixSpecifications
{
    MatrixSpecifications(
        std::size_t const nrows_, std::size_t const ncols_,
        std::vector<GlobalIndexType> const* const ghost_indices_,
        GlobalSparsityPattern const* const sparsity_pattern_)
        : nrows(nrows_),
          ncols(ncols_),
          ghost_indices(ghost_indices_),
          sparsity_pattern(sparsity_pattern_)
    {
    }

    std::size_t const nrows;
    std::size_t const ncols;
    std::vector<GlobalIndexType> const* const ghost_indices;
    GlobalSparsityPattern const* const sparsity_pattern;
};

}  // namespace MathLib
