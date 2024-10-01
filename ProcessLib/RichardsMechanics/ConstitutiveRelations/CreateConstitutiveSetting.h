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
        BaseLib::ConfigTree const& config);
};
}  // namespace ProcessLib::RichardsMechanics
