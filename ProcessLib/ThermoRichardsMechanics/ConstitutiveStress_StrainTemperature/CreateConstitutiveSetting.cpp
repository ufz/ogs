/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateConstitutiveSetting.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
std::map<int, std::shared_ptr<SolidConstitutiveRelation<DisplacementDim>>>
CreateConstitutiveSetting<DisplacementDim>::createSolidConstitutiveRelations(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    MeshLib::PropertyVector<int> const* const material_ids,
    BaseLib::ConfigTree const& config)
{
    return MaterialLib::Solids::createConstitutiveRelations<DisplacementDim>(
        parameters, local_coordinate_system, material_ids, config);
}

template struct CreateConstitutiveSetting<2>;
template struct CreateConstitutiveSetting<3>;
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
