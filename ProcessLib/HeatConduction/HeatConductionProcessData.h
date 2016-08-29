/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_HEATCONDUCTION_HEATCONDUCTIONPROCESSDATA_H
#define PROCESSLIB_HEATCONDUCTION_HEATCONDUCTIONPROCESSDATA_H

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename ReturnType, typename... Args>
struct Parameter;

namespace HeatConduction
{
// copy that and change names
struct HeatConductionProcessData
{
    HeatConductionProcessData(
        ProcessLib::Parameter<double, MeshLib::Element const&> const&
            thermal_conductivity_,
        ProcessLib::Parameter<double, MeshLib::Element const&> const&
            heat_capacity_,
        ProcessLib::Parameter<double, MeshLib::Element const&> const& density_)
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

    Parameter<double, MeshLib::Element const&> const& thermal_conductivity;
    Parameter<double, MeshLib::Element const&> const& heat_capacity;
    Parameter<double, MeshLib::Element const&> const& density;
};

}  // namespace HeatConduction
}  // namespace ProcessLib

#endif  // PROCESSLIB_HEATCONDUCTION_HEATCONDUCTIONPROCESSDATA_H
