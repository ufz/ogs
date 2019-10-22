/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <unordered_map>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MeshLib/PropertyVector.h"
#include "ProcessLib/BoundaryCondition/Python/BHEInflowPythonBoundaryConditionPythonSideInterface.h"
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"
namespace MeshLib
{
class Element;
}

namespace ProcessLib::HeatTransportBHE
{
struct HeatTransportBHEProcessData final
{
    HeatTransportBHEProcessData(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        std::vector<BHE::BHETypes>&& vec_BHEs_,
        bool const& has_network_python_bc_,
        BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object_ =
            nullptr)
        : media_map(std::move(media_map_)),
          _vec_BHE_property(std::move(vec_BHEs_)),
          has_network_python_bc(has_network_python_bc_),
          py_bc_object(py_bc_object_)
    {
    }
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::unordered_map<int, int> _map_materialID_to_BHE_ID{};
    // if bhe network exists python boundary condition
    bool has_network_python_bc;
    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object;
};
}  // namespace ProcessLib::HeatTransportBHE
