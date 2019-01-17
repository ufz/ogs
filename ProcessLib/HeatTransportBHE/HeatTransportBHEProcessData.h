/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <unordered_map>

#include "MeshLib/PropertyVector.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"

#include "ProcessLib/BoundaryCondition/Python/BHEInflowPythonBoundaryConditionPythonSideInterface.h"

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
        std::vector<BHE::BHETypes>&& vec_BHEs_,
        bool const& if_bhe_network_exist_python_bc_,
        BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object_)
        : thermal_conductivity_solid(thermal_conductivity_solid_),
          thermal_conductivity_fluid(thermal_conductivity_fluid_),
          thermal_conductivity_gas(thermal_conductivity_gas_),
          heat_capacity_solid(heat_capacity_solid_),
          heat_capacity_fluid(heat_capacity_fluid_),
          heat_capacity_gas(heat_capacity_gas_),
          density_solid(density_solid_),
          density_fluid(density_fluid_),
          density_gas(density_gas_),
          _vec_BHE_property(std::move(vec_BHEs_)),
          if_bhe_network_exist_python_bc(if_bhe_network_exist_python_bc_),
          py_bc_object(py_bc_object_)
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
    std::unordered_map<int, int> _map_materialID_to_BHE_ID;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    // if bhe network exists python boundary condition
    bool if_bhe_network_exist_python_bc;

    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object;
};
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
