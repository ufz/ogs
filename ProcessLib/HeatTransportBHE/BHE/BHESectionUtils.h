// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "BaseLib/Error.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/// Compute cumulative arc-length distances from the first node.
inline std::vector<double> cumulativeDistances(
    std::vector<MeshLib::Node*> const& sorted_nodes)
{
    std::vector<double> distances;
    distances.reserve(sorted_nodes.size());
    distances.push_back(0.0);

    double cumulative_distance = 0.0;
    for (std::size_t i = 1; i < sorted_nodes.size(); ++i)
    {
        cumulative_distance +=
            std::sqrt(MathLib::sqrDist(*sorted_nodes[i - 1], *sorted_nodes[i]));
        distances.push_back(cumulative_distance);
    }

    return distances;
}

/// Group consecutive sections with the same diameter value into one section.
/// Returns {section_boundaries, section_diameters}.
/// Uses absolute tolerance for floating-point comparison.
inline std::pair<std::vector<double>, std::vector<double>> groupSections(
    std::vector<double> const& distances, std::vector<double> const& diameters)
{
    if (distances.empty() || diameters.empty() ||
        distances.size() != diameters.size())
    {
        OGS_FATAL(
            "Invalid sampled diameter profile. Distances and diameters "
            "must be non-empty and have the same size.");
    }

    if (std::any_of(diameters.begin(), diameters.end(),
                    [](double d) { return d <= 0.0; }))
    {
        OGS_FATAL("All borehole diameters must be positive.");
    }

    std::vector<double> section_boundaries;
    std::vector<double> section_diameters;
    section_boundaries.push_back(distances[0]);
    section_diameters.push_back(diameters[0]);

    for (std::size_t i = 1; i < diameters.size(); ++i)
    {
        constexpr double grouping_tolerance = 1e-12;
        if (std::abs(diameters[i] - section_diameters.back()) >
            grouping_tolerance)
        {
            section_boundaries.push_back(distances[i]);
            section_diameters.push_back(diameters[i]);
        }
    }

    return {section_boundaries, section_diameters};
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
