/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "BaseLib/StringTools.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib::Solids
{
template <typename SolidConstitutiveRelation>
std::map<int, std::shared_ptr<SolidConstitutiveRelation>>
createConstitutiveRelationsGeneric(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config,
    std::unique_ptr<SolidConstitutiveRelation> (*create_constitutive_relation)(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&,
        std::optional<ParameterLib::CoordinateSystem> const&,
        BaseLib::ConfigTree const&))
{
    auto const constitutive_relation_configs =
        //! \ogs_file_param{material__solid__constitutive_relation}
        config.getConfigSubtreeList("constitutive_relation");

    std::map<int, std::shared_ptr<SolidConstitutiveRelation>>
        constitutive_relations;

    for (auto const& constitutive_relation_config :
         constitutive_relation_configs)
    {
        auto const material_id_string =
            //! \ogs_file_attr{material__solid__constitutive_relation__id}
            constitutive_relation_config.getConfigAttribute<std::string>("id",
                                                                         "0");

        auto const material_ids_of_this_constitutive_relation =
            BaseLib::splitMaterialIdString(material_id_string);

        auto first_relation_for_material_id = constitutive_relations.end();
        for (auto const& material_id :
             material_ids_of_this_constitutive_relation)
        {
            if (constitutive_relations.find(material_id) !=
                constitutive_relations.end())
            {
                OGS_FATAL(
                    "Multiple constitutive relations were specified for the "
                    "same material id {:d}. Keep in mind, that if no material "
                    "id is specified, it is assumed to be 0 by default.",
                    material_id);
            }
            if (material_id == material_ids_of_this_constitutive_relation[0])
            {
                auto [it, insertion_succeeded] = constitutive_relations.emplace(
                    material_id,
                    create_constitutive_relation(parameters,
                                                 local_coordinate_system,
                                                 constitutive_relation_config));
                assert(insertion_succeeded);
                first_relation_for_material_id = it;
            }
            else
            {
                // This constitutive relation has multiple material IDs assigned
                // and this is not the first material ID. Therefore we can reuse
                // the constitutive relation we created before.
                assert(first_relation_for_material_id !=
                       constitutive_relations.end());
                constitutive_relations.emplace(
                    material_id, first_relation_for_material_id->second);
            }
        }
    }

    DBUG("Found {:d} constitutive relations.", constitutive_relations.size());

    return constitutive_relations;
}

}  // namespace MaterialLib::Solids
