/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <spdlog/fmt/bundled/ranges.h>

#include <map>
#include <range/v3/view/map.hpp>

#include "MeshLib/PropertyVector.h"

namespace MaterialLib
{
namespace Solids
{
/// Choose solid material model for given element id out of a set of models,
/// possibly using the material ids.
///
/// Only two possibilities yield a valid result and result in OGS_FATAL call
/// otherwise.
/// 1. If the material id is not defined then search for a constitutive
/// relation with id 0 (the default value if not specified).
/// 2. There is only one constitutive relation with id 0 then the material and
/// element ids are ignored and the only constitutive relation (under id 0) is
/// returned.
/// 3. If material ids are defined then search for a constitutive relation
/// corresponding to the material id of the current element.
template <typename SolidMaterialsMap>
auto& selectSolidConstitutiveRelation(
    SolidMaterialsMap const& constitutive_relations,
    MeshLib::PropertyVector<int> const* const material_ids,
    std::size_t const element_id)
{
    // Multiple constitutive relations and no material ids should not be valid.
    if (constitutive_relations.size() > 1 && material_ids == nullptr)
    {
        OGS_FATAL(
            "There are {} constitutive relations provided in the project file "
            "but no MaterialIDs could be found in the mesh.",
            constitutive_relations.size());
    }

    int const material_id = ((constitutive_relations.size() == 1 &&
                              constitutive_relations.begin()->first == 0) ||
                             material_ids == nullptr)
                                ? 0
                                : (*material_ids)[element_id];

    auto const constitutive_relation = constitutive_relations.find(material_id);
    if (constitutive_relation == end(constitutive_relations))
    {
        OGS_FATAL(
            "No constitutive relation found for material id {:d} and element "
            "{:d}. There are {:d} constitutive relations available, "
            "corresponding to the ids: {}",
            material_id, element_id, constitutive_relations.size(),
            fmt::join(constitutive_relations | ranges::views::keys, " "));
    }

    if (constitutive_relation->second == nullptr)
    {
        OGS_FATAL(
            "The constitutive relation found for material id {:d} and element "
            "{:d} is a nullptr, which is impossible.",
            material_id, element_id);
    }

    return *constitutive_relation->second;
}
}  // namespace Solids
}  // namespace MaterialLib
