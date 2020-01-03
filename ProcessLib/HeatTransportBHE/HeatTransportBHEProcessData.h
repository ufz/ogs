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
#include "ProcessLib/HeatTransportBHE/BHE/BHETypes.h"

namespace MeshLib
{
class Element;
}

namespace ProcessLib::HeatTransportBHE
{
struct HeatTransportBHEProcessData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    std::vector<BHE::BHETypes> _vec_BHE_property;

    MeshLib::PropertyVector<int> const* _mesh_prop_materialIDs = nullptr;
    std::unordered_map<int, int> _map_materialID_to_BHE_ID{};
};
}  // namespace ProcessLib::HeatTransportBHE
