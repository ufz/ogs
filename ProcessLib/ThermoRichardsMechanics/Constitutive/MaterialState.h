/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct MaterialStateData
{
    explicit MaterialStateData(
        MaterialLib::Solids::MechanicsBase<DisplacementDim> const&
            solid_material)
        : material_state_variables(
              solid_material.createMaterialStateVariables())
    {
    }

    void pushBackState() { material_state_variables->pushBackState(); }

    std::unique_ptr<typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables>
        material_state_variables;
};
}  // namespace ProcessLib::ThermoRichardsMechanics
