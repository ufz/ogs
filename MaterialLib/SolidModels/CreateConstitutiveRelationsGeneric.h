/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/Utils/MediaCreation.h"
#include "MeshLib/PropertyVector.h"
#include "ParameterLib/Parameter.h"

namespace MaterialLib::Solids
{
template <typename SolidConstitutiveRelation>
std::map<int, std::shared_ptr<SolidConstitutiveRelation>>
createConstitutiveRelationsGeneric(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    MeshLib::PropertyVector<int> const* const material_ids,
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
        auto create = [&create_constitutive_relation,
                       &parameters,
                       &local_coordinate_system,
                       &constitutive_relation_config](int const /*id*/)
        {
            return create_constitutive_relation(parameters,
                                                local_coordinate_system,
                                                constitutive_relation_config);
        };

        auto const material_id_string =
            //! \ogs_file_attr{material__solid__constitutive_relation__id}
            constitutive_relation_config.getConfigAttribute<std::string>("id",
                                                                         "0");

        auto const material_ids_of_this_constitutive_relation =
            MaterialLib::parseMaterialIdString(material_id_string,
                                               material_ids);

        for (auto const& material_id :
             material_ids_of_this_constitutive_relation)
        {
            MaterialLib::createMediumForId(
                material_id,
                constitutive_relations,
                material_ids_of_this_constitutive_relation,
                create);
        }
    }

    DBUG("Found {:d} constitutive relations.", constitutive_relations.size());

    return constitutive_relations;
}

}  // namespace MaterialLib::Solids
