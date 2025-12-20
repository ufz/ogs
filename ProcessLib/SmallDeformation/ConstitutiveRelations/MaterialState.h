// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::SmallDeformation
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
}  // namespace ProcessLib::SmallDeformation
