// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ParameterLib/Parameter.h"

namespace ProcessLib::RichardsMechanics
{
template <int DisplacementDim>
struct CreateConstitutiveSetting
{
    static std::map<
        int,
        std::shared_ptr<MaterialLib::Solids::MechanicsBase<DisplacementDim>>>
    createSolidConstitutiveRelations(
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        std::optional<ParameterLib::CoordinateSystem> const&
            local_coordinate_system,
        MeshLib::PropertyVector<int> const* const material_ids,
        BaseLib::ConfigTree const& config);
};
}  // namespace ProcessLib::RichardsMechanics
