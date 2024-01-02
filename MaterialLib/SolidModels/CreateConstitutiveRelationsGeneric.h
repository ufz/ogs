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
#include "ParameterLib/Parameter.h"

namespace MaterialLib::Solids
{
template <typename SolidConstitutiveRelation>
std::map<int, std::unique_ptr<SolidConstitutiveRelation>>
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

    std::map<int, std::unique_ptr<SolidConstitutiveRelation>>
        constitutive_relations;

    for (auto const& constitutive_relation_config :
         constitutive_relation_configs)
    {
        int const material_id =
            //! \ogs_file_attr{material__solid__constitutive_relation__id}
            constitutive_relation_config.getConfigAttribute<int>("id", 0);

        if (constitutive_relations.find(material_id) !=
            constitutive_relations.end())
        {
            OGS_FATAL(
                "Multiple constitutive relations were specified for the same "
                "material id {:d}. Keep in mind, that if no material id is "
                "specified, it is assumed to be 0 by default.",
                material_id);
        }

        constitutive_relations.emplace(
            material_id,
            create_constitutive_relation(parameters,
                                         local_coordinate_system,
                                         constitutive_relation_config));
    }

    DBUG("Found {:d} constitutive relations.", constitutive_relations.size());

    return constitutive_relations;
}

}  // namespace MaterialLib::Solids
