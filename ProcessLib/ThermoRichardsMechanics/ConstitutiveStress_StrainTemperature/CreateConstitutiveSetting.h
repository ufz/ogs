// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>

#include "ParameterLib/Parameter.h"
#include "TraitsBase.h"

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
struct CreateConstitutiveSetting
{
    static std::map<int,
                    std::shared_ptr<SolidConstitutiveRelation<DisplacementDim>>>
    createSolidConstitutiveRelations(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system,
        MeshLib::PropertyVector<int> const* const material_ids,
        BaseLib::ConfigTree const& config);
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
