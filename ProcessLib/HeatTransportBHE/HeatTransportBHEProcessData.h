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
        BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object_ =
            nullptr)
        : media_map(std::move(media_map_)),
          vec_BHE_property_(std::move(vec_BHEs_)),
          py_bc_object(py_bc_object_)
    {
    }
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    std::vector<BHE::BHETypes> vec_BHE_property_;

    MeshLib::PropertyVector<int> const* mesh_prop_materialIDs_ = nullptr;
    std::unordered_map<int, int> map_materialID_to_BHE_ID_{};

    //! Python object computing BC values.
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object;
};
}  // namespace ProcessLib::HeatTransportBHE
