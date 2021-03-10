/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on July 10, 2018, 12:09 PM
 */

#include "CreateConstitutiveRelation.h"

#include "CreateCreepBGRa.h"
#include "CreateEhlers.h"
#include "CreateLinearElasticIsotropic.h"
#include "CreateLinearElasticOrthotropic.h"
#include "CreateLubby2.h"
#ifdef OGS_USE_MFRONT
#include "MFront/CreateMFront.h"
#endif  // OGS_USE_MFRONT

#include "MechanicsBase.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ParameterLib/Parameter.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelation(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto const type =
        //! \ogs_file_param{material__solid__constitutive_relation__type}
        config.peekConfigParameter<std::string>("type");

    if (type == "Ehlers")
    {
        return MaterialLib::Solids::Ehlers::createEhlers<DisplacementDim>(
            parameters, config);
    }
    if (type == "LinearElasticIsotropic")
    {
        const bool skip_type_checking = false;
        return MaterialLib::Solids::createLinearElasticIsotropic<
            DisplacementDim>(parameters, config, skip_type_checking);
    }
    if (type == "LinearElasticOrthotropic")
    {
        const bool skip_type_checking = false;
        return MaterialLib::Solids::createLinearElasticOrthotropic<
            DisplacementDim>(
            parameters, local_coordinate_system, config, skip_type_checking);
    }
    if (type == "Lubby2")
    {
        return MaterialLib::Solids::Lubby2::createLubby2<DisplacementDim>(
            parameters, config);
    }
    if (type == "CreepBGRa")
    {
        return MaterialLib::Solids::Creep::createCreepBGRa<DisplacementDim>(
            parameters, config);
    }
    if (type == "MFront")
    {
#ifdef OGS_USE_MFRONT
        return MaterialLib::Solids::MFront::createMFront<DisplacementDim>(
            parameters, local_coordinate_system, config);
#else   // OGS_USE_MFRONT
        OGS_FATAL(
            "OGS is compiled without MFront support. See OGS_USE_MFRONT CMake "
            "option.");
#endif  // OGS_USE_MFRONT
    }
    OGS_FATAL("Cannot construct constitutive relation of given type '{:s}'.",
              type);
}

template <int DisplacementDim>
std::map<int,
         std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
createConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto const constitutive_relation_configs =
        //! \ogs_file_param{material__solid__constitutive_relation}
        config.getConfigSubtreeList("constitutive_relation");

    std::map<int, std::unique_ptr<MechanicsBase<DisplacementDim>>>
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
            createConstitutiveRelation<DisplacementDim>(
                parameters,
                local_coordinate_system,
                constitutive_relation_config));
    }

    DBUG("Found {:d} constitutive relations.", constitutive_relations.size());

    return constitutive_relations;
}

template std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>>
createConstitutiveRelations<2>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template std::map<int, std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>>
createConstitutiveRelations<3>(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace Solids
}  // namespace MaterialLib
