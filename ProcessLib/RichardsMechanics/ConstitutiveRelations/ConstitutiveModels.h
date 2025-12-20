// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"
#include "ProcessLib/Graph/ConstructModels.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/SpecificBodyForceData.h"

namespace ProcessLib::RichardsMechanics
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
using ConstitutiveModels = std::tuple<>;

template <int DisplacementDim, typename TRMProcessData>
ConstitutiveModels<DisplacementDim> createConstitutiveModels(
    TRMProcessData const& process_data,
    MaterialLib::Solids::MechanicsBase<DisplacementDim> const& solid_material)
{
    return ProcessLib::Graph::constructModels<
        ConstitutiveModels<DisplacementDim>>(
        ProcessLib::ThermoRichardsMechanics::SpecificBodyForceData<
            DisplacementDim>{process_data.specific_body_force},
        solid_material);
}
}  // namespace ProcessLib::RichardsMechanics
