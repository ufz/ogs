/**
 * \file
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
#include "ParameterLib/Parameter.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace MeshLib
{
class Element;
}

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib::HeatTransportBHE
{
struct HeatTransportBHEProcessData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    // ! thermal conductivity values for the three phases
    ParameterLib::Parameter<double> const& thermal_conductivity_solid;
    ParameterLib::Parameter<double> const& thermal_conductivity_fluid;
    ParameterLib::Parameter<double> const& thermal_conductivity_gas;

    // ! heat capacity values for the three phases
    ParameterLib::Parameter<double> const& heat_capacity_solid;
    ParameterLib::Parameter<double> const& heat_capacity_fluid;
    ParameterLib::Parameter<double> const& heat_capacity_gas;

    // ! density values for the three phases
    ParameterLib::Parameter<double> const& density_solid;
    ParameterLib::Parameter<double> const& density_fluid;
    ParameterLib::Parameter<double> const& density_gas;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::unordered_map<int, int> _map_materialID_to_BHE_ID{};
};
}  // namespace ProcessLib::HeatTransportBHE
