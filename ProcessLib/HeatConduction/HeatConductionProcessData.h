/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib::HeatConduction
{
struct HeatConductionProcessData
{
    HeatConductionProcessData(
        ParameterLib::Parameter<double> const& thermal_conductivity_,
        ParameterLib::Parameter<double> const& heat_capacity_,
        ParameterLib::Parameter<double> const& density_)
        : thermal_conductivity(thermal_conductivity_),
          heat_capacity(heat_capacity_),
          density(density_)
    {
    }

    HeatConductionProcessData(HeatConductionProcessData&& other)
        : thermal_conductivity(other.thermal_conductivity),
          heat_capacity(other.heat_capacity),
          density(other.density)
    {
    }

    //! Copies are forbidden.
    HeatConductionProcessData(HeatConductionProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatConductionProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatConductionProcessData&&) = delete;

    ParameterLib::Parameter<double> const& thermal_conductivity;
    ParameterLib::Parameter<double> const& heat_capacity;
    ParameterLib::Parameter<double> const& density;
};
}  // namespace ProcessLib::HeatConduction
