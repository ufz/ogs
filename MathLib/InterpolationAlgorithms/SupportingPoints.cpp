// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SupportingPoints.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <range/v3/algorithm/find_if.hpp>

#include "BaseLib/Error.h"
#include "BaseLib/quicksort.h"

namespace MathLib
{
void sortAndCheckSupportingPoints(std::vector<double>& supporting_points,
                                  std::vector<double>& values_at_supp_pnts,
                                  std::size_t const min_size,
                                  std::string_view const description,
                                  bool const supp_pnts_sorted)
{
    if (supporting_points.size() != values_at_supp_pnts.size())
    {
        OGS_FATAL("The {:s} has {:d} supporting points but {:d} values.",
                  description, supporting_points.size(),
                  values_at_supp_pnts.size());
    }
    if (supporting_points.size() < min_size)
    {
        OGS_FATAL(
            "The {:s} needs at least {:d} supporting points, but only {:d} "
            "were given.",
            description, min_size, supporting_points.size());
    }
    // NaN compares false with everything, so it would pass the sorting and
    // the check below unnoticed and leave the points unsorted.
    if (auto const it =
            ranges::find_if(supporting_points,
                            [](double const x) { return !std::isfinite(x); });
        it != supporting_points.end())
    {
        OGS_FATAL("The supporting point {:d} ({:g}) of the {:s} is not finite.",
                  std::distance(supporting_points.begin(), it), *it,
                  description);
    }
    if (!supp_pnts_sorted)
    {
        BaseLib::quicksort(supporting_points.begin(), supporting_points.end(),
                           values_at_supp_pnts.begin());
    }
    // Finds equal supporting points and, if sorting was skipped, also
    // decreasing ones.
    if (auto const it = std::adjacent_find(supporting_points.begin(),
                                           supporting_points.end(),
                                           std::greater_equal<>{});
        it != supporting_points.end())
    {
        auto const i = std::distance(supporting_points.begin(), it);
        OGS_FATAL(
            "The supporting points {:d} ({:g}) and {:d} ({:g}) of the {:s} are "
            "not strictly increasing.",
            i, *it, i + 1, *std::next(it), description);
    }
}
}  // namespace MathLib
