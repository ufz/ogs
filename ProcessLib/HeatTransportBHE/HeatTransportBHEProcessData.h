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
#ifdef OGS_USE_PYTHON
    HeatTransportBHEProcessData(
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        std::vector<BHE::BHETypes>&& vec_BHEs_,
        bool const& if_bhe_network_exist_python_bc_,
        BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object_)
        : media_map(std::move(media_map_)),
          _vec_BHE_property(std::move(vec_BHEs_)),
          if_bhe_network_exist_python_bc(if_bhe_network_exist_python_bc_),
          py_bc_object(py_bc_object_)
#else
        std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>&&
            media_map_,
        std::vector<BHE::BHETypes>&& vec_BHEs_,
        bool const& if_bhe_network_exist_python_bc_)
        : media_map(std::move(media_map_)),
          _vec_BHE_property(std::move(vec_BHEs_)),
          if_bhe_network_exist_python_bc(if_bhe_network_exist_python_bc_)
#endif  // OGS_USE_PYTHON
    {
    }
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::unordered_map<int, int> _map_materialID_to_BHE_ID{};
    // if bhe network exists python boundary condition
    bool if_bhe_network_exist_python_bc;
#ifdef OGS_USE_PYTHON
    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object;
#endif
};
}  // namespace ProcessLib::HeatTransportBHE
