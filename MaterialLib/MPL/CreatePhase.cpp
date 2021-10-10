/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{
    using namespace MaterialPropertyLib;

    //! \ogs_file_param{prj__media__medium__phases__phase__type}
    auto&& phase_type = config.getConfigParameter<std::string>("type");

    if (phase_type.empty())
    {
        OGS_FATAL("Phase type is a mandatory field and cannot be empty.");
    }

    std::array<std::string, 4> const allowed_phase_types = {
        {"Solid", "AqueousLiquid", "NonAqueousLiquid", "Gas"}};

    if (std::none_of(allowed_phase_types.begin(),
                     allowed_phase_types.end(),
                     [&phase_type](std::string const& type)
                     { return phase_type == type; }))
    {
        ERR("Phase type should be one of:");
        for (auto const& type : allowed_phase_types)
        {
            ERR("{:s}", type);
        }
        OGS_FATAL("Wrong phase type '{:s}' given.", phase_type);
    }

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
            phase_type);
    }

    return std::make_unique<Phase>(std::move(phase_type), std::move(components),
                                   std::move(properties));
}
}  // namespace

namespace MaterialPropertyLib
{
std::vector<std::unique_ptr<Phase>> createPhases(
    int const geometry_dimension,
    std::optional<BaseLib::ConfigTree> const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
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
                         [phase_name = phase->name](auto const& p)
                         { return p->name == phase_name; }) != phases.end())
        {
            OGS_FATAL("Found duplicates with the same phase name tag '{:s}'.",
                      phase->name);
        }

        phases.push_back(std::move(phase));
    }

    return phases;
}
}  // namespace MaterialPropertyLib
