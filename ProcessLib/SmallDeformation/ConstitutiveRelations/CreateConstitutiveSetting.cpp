// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateConstitutiveSetting.h"

#include "MaterialLib/SolidModels/CreateConstitutiveRelation.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
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
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
