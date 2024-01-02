/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateConstitutiveSetting.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
std::map<int, std::unique_ptr<SolidConstitutiveRelation<DisplacementDim>>>
CreateConstitutiveSetting<DisplacementDim>::createSolidConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    return MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
        parameters, local_coordinate_system, config);
}

template struct CreateConstitutiveSetting<2>;
template struct CreateConstitutiveSetting<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
