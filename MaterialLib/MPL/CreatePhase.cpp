// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreatePhase.h"

#include <set>
#include <string>

#include "BaseLib/ConfigTree.h"
#include "CreateComponent.h"
#include "CreateProperty.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "ParameterLib/Parameter.h"
#include "Phase.h"

namespace
{
std::unique_ptr<MaterialPropertyLib::Phase> createPhase(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    using namespace MaterialPropertyLib;

    //! \ogs_file_param{prj__media__medium__phases__phase__type}
    auto&& phase_type_string = config.getConfigParameter<std::string>("type");

    if (phase_type_string.empty())
    {
        OGS_FATAL("Phase type is a mandatory field and cannot be empty.");
    }

    //! \ogs_file_param_special{prj__media__medium__phases__phase__Solid}
    //! \ogs_file_param_special{prj__media__medium__phases__phase__FrozenLiquid}
    //! \ogs_file_param_special{prj__media__medium__phases__phase__AqueousLiquid}
    //! \ogs_file_param_special{prj__media__medium__phases__phase__Gas}
    auto const phase_type = fromString(phase_type_string);

    // Parsing of optional components.
    auto components = createComponents(
        geometry_dimension,
        //! \ogs_file_param{prj__media__medium__phases__phase__components}
        config.getConfigSubtreeOptional("components"), parameters,
        local_coordinate_system, curves);

    // Properties of optional properties.
    auto properties = createProperties(
        geometry_dimension,
        //! \ogs_file_param{prj__media__medium__phases__phase__properties}
        config.getConfigSubtreeOptional("properties"), parameters,
        local_coordinate_system, curves);

    if (components.empty() && !properties)
    {
        OGS_FATAL(
            "Neither tag <components> nor tag <properties> has been set for "
            "the phase '{:s}'.",
            toString(phase_type));
    }

    return std::make_unique<Phase>(phase_type, std::move(components),
                                   std::move(properties));
}
}  // namespace

namespace MaterialPropertyLib
{
std::vector<std::unique_ptr<Phase>> createPhases(
    int const geometry_dimension,
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    if (!config)
    {
        return {};
    }

    std::vector<std::unique_ptr<Phase>> phases;

    for (auto phase_config :
         //! \ogs_file_param{prj__media__medium__phases__phase}
         config->getConfigSubtreeList("phase"))
    {
        auto phase = createPhase(geometry_dimension, phase_config, parameters,
                                 local_coordinate_system, curves);

        if (std::find_if(phases.begin(),
                         phases.end(),
                         [phase_type = phase->phaseName](auto const& p)
                         { return p->phaseName == phase_type; }) !=
            phases.end())
        {
            OGS_FATAL("Found duplicates with the same phase name tag '{:s}'.",
                      toString(phase->phaseName));
        }

        phases.push_back(std::move(phase));
    }

    return phases;
}
}  // namespace MaterialPropertyLib
