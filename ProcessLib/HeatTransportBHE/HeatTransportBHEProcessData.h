/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/PropertyVector.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib
{
template <typename T>
struct Parameter;

namespace HeatTransportBHE
{
using namespace BHE;

struct HeatTransportBHEProcessData
{
    HeatTransportBHEProcessData(
        Parameter<double> const& thermal_conductivity_solid_,
        Parameter<double> const& thermal_conductivity_fluid_,
        Parameter<double> const& thermal_conductivity_gas_,
        Parameter<double> const& heat_capacity_solid_,
        Parameter<double> const& heat_capacity_fluid_,
        Parameter<double> const& heat_capacity_gas_,
        Parameter<double> const& density_solid_,
        Parameter<double> const& density_fluid_,
        Parameter<double> const& density_gas_,
        std::vector<std::unique_ptr<BHE::BHEAbstract>>&& vec_BHEs_)
        : thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          thermal_conductivity_gas(thermal_conductivity_gas_),
          heat_capacity_solid(heat_capacity_solid_),
          heat_capacity_fluid(heat_capacity_fluid_),
          heat_capacity_gas(heat_capacity_gas_),
          density_solid(density_solid_),
          density_fluid(density_fluid_),
          density_gas(density_gas_),
          _vec_BHE_property(std::move(vec_BHEs_))
    {
    }

    HeatTransportBHEProcessData(HeatTransportBHEProcessData&& other) = default;

    //! Copies are forbidden.
    HeatTransportBHEProcessData(HeatTransportBHEProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatTransportBHEProcessData const&) = delete;

    //! Assignments are not needed.
    void operator=(HeatTransportBHEProcessData&&) = delete;

    // ! thermal conductivity values for the three phases
    Parameter<double> const& thermal_conductivity_solid;
    Parameter<double> const& thermal_conductivity_fluid;
    Parameter<double> const& thermal_conductivity_gas;

    // ! heat capacity values for the three phases
    Parameter<double> const& heat_capacity_solid;
    Parameter<double> const& heat_capacity_fluid;
    Parameter<double> const& heat_capacity_gas;

    // ! density values for the three phases
    Parameter<double> const& density_solid;
    Parameter<double> const& density_fluid;
    Parameter<double> const& density_gas;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::vector<std::size_t> _map_materialID_to_BHE_ID;

    std::vector<std::unique_ptr<BHE::BHEAbstract>> _vec_BHE_property;

    // a table of connected BHE IDs for each element
    std::vector<std::vector<std::size_t>> _vec_ele_connected_BHE_IDs;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
