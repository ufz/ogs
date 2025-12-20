// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace MathLib
{
/// A vector telling how many nonzeros there are in each global matrix row.
template <typename IndexType>
using SparsityPattern = std::vector<IndexType>;
}  // namespace MathLib
