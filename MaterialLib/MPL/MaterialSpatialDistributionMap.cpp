/**
 * \file
 * \date   Nov 28, 2017
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "MaterialSpatialDistributionMap.h"

#include <spdlog/fmt/ranges.h>

#include <range/v3/view/map.hpp>

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

    assert(!media_.empty());

    if (auto const it = media_.find(material_id); it != media_.end())
    {
        return it->second.get();
    }

    //
    // Error handling until end of the function.
    //

    if (material_ids_ == nullptr)
    {
        assert(material_id == 0);
        ERR("No MaterialIDs given in the mesh therefore default material id = "
            "0 is used.");
    }
    auto keys = media_ | ranges::views::keys;

    if (media_.size() == 1)
    {
        ERR("Single medium for material id {} is defined.",
            fmt::join(keys, ", "));
    }
    else
    {
        ERR("Media for material ids {} are defined.", fmt::join(keys, ", "));
    }
    OGS_FATAL("No medium for material id {} found for element {}.", material_id,
              element_id);
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
