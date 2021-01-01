/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>

#include "MeshLib/PropertyVector.h"

#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
/// Choose solid material model for given element id out of a set of models,
/// possibly using the material ids.
///
/// Only two possibilities yield a valid result and result in OGS_FATAL call
/// otherwise.
/// 1. There is only one constitutive relation in the set of
/// constitutive_relations, then the material and element ids are ignored and
/// the only constitutive relation (under id 0) is returned.
/// 2. For each material id chosen for the given element id there exists a
/// constitutive relation in the set.
template <int DisplacementDim>
MechanicsBase<DisplacementDim>& selectSolidConstitutiveRelation(
    std::map<int, std::unique_ptr<MechanicsBase<DisplacementDim>>> const&
        constitutive_relations,
    MeshLib::PropertyVector<int> const* const material_ids,
    std::size_t const element_id)
{
    int material_id;
    if (constitutive_relations.size() == 1 || material_ids == nullptr)
    {
        material_id = 0;
    }
    else
    {
        material_id = (*material_ids)[element_id];
    }

    auto const constitutive_relation = constitutive_relations.find(material_id);
    if (constitutive_relation == end(constitutive_relations))
    {
        OGS_FATAL(
            "No constitutive relation found for material id {:d} and element "
            "{:d}. There are {:d} constitutive relations available.",
            material_id, element_id, constitutive_relations.size());
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
