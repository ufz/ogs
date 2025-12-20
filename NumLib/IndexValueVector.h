// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

namespace NumLib
{
template <typename IndexType>
struct IndexValueVector final
{
    std::vector<IndexType> ids;
    std::vector<double> values;
};

}  // namespace NumLib
