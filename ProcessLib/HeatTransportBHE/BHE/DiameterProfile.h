// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <algorithm>
#include <functional>
#include <numbers>
#include <vector>

#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct DiameterProfile
{
    DiameterProfile(std::vector<double> bounds, std::vector<double> diams)
        : section_boundaries(std::move(bounds)),
          section_diameters(std::move(diams))
    {
        if (section_boundaries.empty() || section_diameters.empty())
        {
            OGS_FATAL("DiameterProfile: vectors must be non-empty.");
        }
        if (section_boundaries.size() != section_diameters.size())
        {
            OGS_FATAL(
                "DiameterProfile: section_boundaries (size {:d}) and "
                "section_diameters (size {:d}) must have the same size.",
                section_boundaries.size(), section_diameters.size());
        }
        if (std::any_of(section_diameters.begin(), section_diameters.end(),
                        [](double d) { return d <= 0.0; }))
        {
            OGS_FATAL("DiameterProfile: all diameters must be positive.");
        }
        if (std::adjacent_find(section_boundaries.begin(),
                               section_boundaries.end(),
                               std::greater_equal<>{}) !=
            section_boundaries.end())
        {
            OGS_FATAL(
                "DiameterProfile: section_boundaries must be strictly "
                "increasing.");
        }
    }

    /// Cumulative distances from wellhead [m] marking section boundaries.
    /// section_boundaries[i] marks the start of section i.
    /// Length = number of sections.
    std::vector<double> const section_boundaries;

    /// Diameters [m] for each section.
    /// section_diameters[i] is the diameter for section i.
    /// Length = number of sections.
    std::vector<double> const section_diameters;

    /// \brief Find the section index for a given distance from wellhead.
    /// Uses binary search for O(log n) lookup.
    /// \param distance_from_wellhead Distance from wellhead [m].
    /// \return Section index (0-based). Returns last section if distance
    /// exceeds last boundary.
    int getSectionIndex(double distance_from_wellhead) const;

    /// \brief Get diameter at a given distance from wellhead.
    /// \param distance_from_wellhead Distance from wellhead [m].
    /// \return Diameter [m] for the corresponding section.
    double diameterAtDistance(double distance_from_wellhead) const
    {
        return section_diameters[getSectionIndex(distance_from_wellhead)];
    }

    /// \brief Get cross-sectional area at a given distance from wellhead.
    /// \param distance_from_wellhead Distance from wellhead [m].
    /// \return Cross-sectional area [m^2] for the corresponding section.
    double areaAtDistance(double distance_from_wellhead) const
    {
        return circleArea(diameterAtDistance(distance_from_wellhead));
    }

    /// \brief Get diameter for a specific section.
    /// \param section_index Section index (0-based).
    /// \return Diameter [m].
    double diameterAtSection(int const section_index) const
    {
        if (section_index < 0 ||
            section_index >= static_cast<int>(section_diameters.size()))
        {
            OGS_FATAL("Invalid section index: {:d}", section_index);
        }
        return section_diameters[section_index];
    }

    /// \brief Get cross-sectional area for a specific section.
    /// \param section_index Section index (0-based).
    /// \return Cross-sectional area [m²].
    double areaAtSection(int const section_index) const
    {
        return circleArea(diameterAtSection(section_index));
    }

    /// \brief Get the number of sections.
    /// \return Number of sections.
    int getNumberOfSections() const
    {
        return static_cast<int>(section_diameters.size());
    }

private:
    static double circleArea(double const diameter)
    {
        return std::numbers::pi * diameter * diameter / 4;
    }
};

inline int DiameterProfile::getSectionIndex(double distance_from_wellhead) const
{
    if (distance_from_wellhead < section_boundaries.front())
    {
        OGS_FATAL(
            "Distance from wellhead is negative or before first section: "
            "{:g}",
            distance_from_wellhead);
    }

    // Binary search to find the section.
    // Find the first boundary that is greater than distance_from_wellhead.
    auto it = std::upper_bound(section_boundaries.begin(),
                               section_boundaries.end(),
                               distance_from_wellhead);

    // Move back one position to get the section containing this distance.
    if (it != section_boundaries.begin())
    {
        --it;
    }

    return static_cast<int>(std::distance(section_boundaries.begin(), it));
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
