// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

/**
 * \file
 *
 * 1) Diersch_2011_CG
 * Two very important references to understand this class implementations are:
 * Diersch, H-JG, D. Bauer, W. Heidemann, Wolfram Rühaak, and Peter Schätzl.
 * Finite element modeling of borehole heat exchanger systems:
 * Part 1. Fundamentals, Computers & Geosciences,
 * Volume 37, Issue 8, August 2011, Pages 1122-1135, ISSN 0098-3004,
 * http://dx.doi.org/10.1016/j.cageo.2010.08.003.
 *
 * 2) FEFLOW_2014_Springer
 * FEFLOW: Finite Element Modeling of Flow, Mass and Heat Transport in Porous
 * and Fractured Media Diersch, Hans-Joerg, 2014, XXXV, 996 p, Springer.
 */

#pragma once

#include <algorithm>
#include <vector>

#include "BaseLib/Error.h"
#include "BoreholeGeometry.h"
#include "FlowAndTemperatureControl.h"
#include "GroutParameters.h"
#include "RefrigerantProperties.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/// Validate that a grout cross-sectional area (borehole area fraction minus
/// pipe outside area) is positive.  Returns the area on success, calls
/// OGS_FATAL on failure.
inline double checkedGroutArea(double const borehole_area_fraction,
                               double const pipe_outside_area,
                               int const section_index)
{
    double const grout_area = borehole_area_fraction - pipe_outside_area;
    if (grout_area <= 0)
    {
        OGS_FATAL(
            "Non-positive grout cross-sectional area at section {:d}. "
            "Borehole diameter is too small for the pipe dimensions.",
            section_index);
    }
    return grout_area;
}

class BHECommon
{
public:
    BHECommon(BoreholeGeometry const& borehole_geometry_,
              RefrigerantProperties const& refrigerant_,
              GroutParameters const& grout_,
              FlowAndTemperatureControl const& flowAndTemperatureControl_,
              bool const use_python_bcs_)
        : borehole_geometry(borehole_geometry_),
          refrigerant(refrigerant_),
          grout(grout_),
          flowAndTemperatureControl(flowAndTemperatureControl_),
          use_python_bcs(use_python_bcs_)
    {
    }

    BoreholeGeometry const borehole_geometry;
    RefrigerantProperties const refrigerant;
    GroutParameters const grout;
    FlowAndTemperatureControl const flowAndTemperatureControl;
    bool const use_python_bcs;
    constexpr bool isPowerBC() const
    {
        return std::visit([](auto const& ftc) { return ftc.is_power_bc; },
                          flowAndTemperatureControl);
    }

    /// \brief Get number of sections in the borehole geometry.
    int getNumberOfSections() const
    {
        return borehole_geometry.sections.getNumberOfSections();
    }

    /// \brief Get section boundaries (cumulative distance from wellhead).
    std::vector<double> const& getSectionBoundaries() const
    {
        return borehole_geometry.sections.section_boundaries;
    }

    /// \brief Get thermal resistance for a specific section and unknown.
    double thermalResistanceAtSection(int const unknown_index,
                                      int const section_index) const
    {
        auto const& thermal_resistances =
            getThermalResistancesAtSection(section_index);
        if (unknown_index < 0 ||
            unknown_index >= static_cast<int>(thermal_resistances.size()))
        {
            OGS_FATAL("Invalid unknown index {:d} for section {:d}.",
                      unknown_index, section_index);
        }
        return thermal_resistances[unknown_index];
    }

protected:
    /// \brief Recompute _sectional_thermal_resistances for all sections.
    ///
    /// \param calcForSection  Callable \c (int section_index) ->
    /// std::vector<double>
    ///                        that returns the resistance vector for that
    ///                        section.
    template <typename Fn>
    void recomputeSectionalResistances(Fn&& calcForSection)
    {
        int const n = getNumberOfSections();
        _sectional_thermal_resistances.clear();
        _sectional_thermal_resistances.reserve(n);
        for (int i = 0; i < n; ++i)
        {
            _sectional_thermal_resistances.push_back(calcForSection(i));
        }
    }
    /// Flow velocities per section [m/s]. For single-section pipes, will have
    /// one element.
    std::vector<double> _flow_velocities;

    /// \brief Get velocity for a section, clamping to last section if index
    /// exceeds the number of velocity entries. Returns 0 when empty.
    double getClampedFlowVelocity(int const section_index) const
    {
        return _flow_velocities.empty()
                   ? 0.0
                   : _flow_velocities[std::min(
                         section_index,
                         static_cast<int>(_flow_velocities.size()) - 1)];
    }

    /// Thermal resistances per section. Each element is a vector of resistances
    /// for that section's unknowns.
    std::vector<std::vector<double>> _sectional_thermal_resistances;

    /// \brief Get the thermal resistance vector for a given section.
    std::vector<double> const& getThermalResistancesAtSection(
        int const section_index) const
    {
        if (section_index < 0 ||
            section_index >=
                static_cast<int>(_sectional_thermal_resistances.size()))
        {
            OGS_FATAL("Invalid section index: {:d}", section_index);
        }
        return _sectional_thermal_resistances[section_index];
    }
};
}  // end of namespace BHE
}  // end of namespace HeatTransportBHE
}  // end of namespace ProcessLib
