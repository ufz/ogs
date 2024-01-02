/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::LargeDeformation
{
template <int DisplacementDim>
class MaterialStateData
{
    using MSV = typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables;

public:
    explicit MaterialStateData(std::unique_ptr<MSV>&& material_state_variables)
        : material_state_variables(std::move(material_state_variables))
    {
    }

    void pushBackState() { material_state_variables->pushBackState(); }

    std::unique_ptr<MSV> material_state_variables;
};
}  // namespace ProcessLib::LargeDeformation
