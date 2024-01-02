/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateConstitutiveRelationIce.h"

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"
#include "MechanicsBase.h"

namespace MaterialLib
{
namespace Solids
{
template <int DisplacementDim>
std::unique_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto const ice_constitutive_relation_config =
        //! \ogs_file_param{material__solid__ice_constitutive_relation}
        config.getConfigSubtreeOptional("ice_constitutive_relation");

    if (!ice_constitutive_relation_config)
    {
        return nullptr;
    }

    return MaterialLib::Solids::createConstitutiveRelation<DisplacementDim>(
        parameters, local_coordinate_system, *ice_constitutive_relation_config);
}

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<2>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);

template std::unique_ptr<MaterialLib::Solids::MechanicsBase<3>>
createConstitutiveRelationIce(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config);
}  // namespace Solids
}  // namespace MaterialLib
