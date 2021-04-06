/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    return const_cast<Medium*>(
        static_cast<MaterialSpatialDistributionMap const&>(*this).getMedium(
            element_id));
}

Medium const* MaterialSpatialDistributionMap::getMedium(
    std::size_t const element_id) const
{
    auto const material_id =
        material_ids_ == nullptr ? 0 : (*material_ids_)[element_id];

    return media_.at(material_id).get();
}

void MaterialSpatialDistributionMap::checkElementHasMedium(
    std::size_t const element_id) const
{
    auto const material_id =
        material_ids_ == nullptr ? 0 : (*material_ids_)[element_id];
    if (media_.find(material_id) == media_.end())
    {
        OGS_FATAL(
            "There is no medium definition for element {:d} with material "
            "ID {:d}. Please define a medium for each material.",
            element_id, material_id);
    }
}

}  // namespace MaterialPropertyLib
