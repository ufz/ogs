/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "MaterialSpatialDistributionMap.h"
#include "MeshLib/Mesh.h"

namespace MaterialPropertyLib
{
Medium* MaterialSpatialDistributionMap::getMedium(std::size_t const element_id)
{
    auto const material_id =
        _material_ids == nullptr ? 0 : (*_material_ids)[element_id];

    return _media.at(material_id).get();
}
}  // namespace MaterialPropertyLib
