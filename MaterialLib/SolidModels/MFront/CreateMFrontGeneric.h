// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/ConfigTree-fwd.h"
#include "MFrontGeneric.h"

namespace MaterialLib::Solids::MFront
{

struct MFrontConfig
{
    mgis::behaviour::Behaviour behaviour;
    std::vector<ParameterLib::Parameter<double> const*> material_properties;
    std::map<std::string, ParameterLib::Parameter<double> const*>
        state_variables_initial_properties;
};

MFrontConfig createMFrontConfig(
    int const displacement_dim,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    BaseLib::ConfigTree const& config);

template <int DisplacementDim, typename Gradients, typename TDynForces,
          typename ExtStateVars>
std::unique_ptr<
    MFrontGeneric<DisplacementDim, Gradients, TDynForces, ExtStateVars>>
createMFrontGeneric(
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    std::optional<ParameterLib::CoordinateSystem> const&
        local_coordinate_system,
    BaseLib::ConfigTree const& config)
{
    auto conf = createMFrontConfig(DisplacementDim, parameters, config);

    return std::make_unique<
        MFrontGeneric<DisplacementDim, Gradients, TDynForces, ExtStateVars>>(
        std::move(conf.behaviour), std::move(conf.material_properties),
        std::move(conf.state_variables_initial_properties),
        local_coordinate_system);
}
}  // namespace MaterialLib::Solids::MFront
