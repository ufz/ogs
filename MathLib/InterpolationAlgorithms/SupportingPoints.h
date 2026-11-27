// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>
#include <string_view>
#include <vector>

namespace MathLib
{
/// Sorts the supporting points in ascending order, reordering the values
/// accordingly, and checks that there are as many values as supporting
/// points, that there are at least \c min_size of them, and that they are
/// finite and pairwise distinct. Each violation is a fatal error.
/// \param supporting_points Abscissae of the supporting points.
/// \param values_at_supp_pnts Values at the supporting points.
/// \param min_size Minimal number of supporting points.
/// \param description Describes the owner of the supporting points in error
/// messages, e.g. "PCHIP curve 'saturation'".
/// \param supp_pnts_sorted If true, sorting is skipped and the supporting
/// points are only checked to be strictly increasing.
void sortAndCheckSupportingPoints(std::vector<double>& supporting_points,
                                  std::vector<double>& values_at_supp_pnts,
                                  std::size_t min_size,
                                  std::string_view description,
                                  bool supp_pnts_sorted);
}  // namespace MathLib
