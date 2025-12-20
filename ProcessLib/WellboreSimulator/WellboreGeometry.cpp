// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "WellboreGeometry.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
WellboreGeometry createWellboreGeometry(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    double const length =
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__length}
        config.getConfigParameter<double>("length");

    auto const& diameter = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__diameter}
        config.getConfigParameter<std::string>("diameter"),
        parameters,
        1,
        nullptr);

    auto const& casing_thickness = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__casing_thickness}
        config.getConfigParameter<std::string>("casing_thickness"),
        parameters,
        1,
        nullptr);

    auto const& pipe_thickness = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__pipe_thickness}
        config.getConfigParameter<std::string>("pipe_thickness"),
        parameters,
        1,
        nullptr);

    auto const& roughness = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__wellbore__roughness}
        config.getConfigParameter<std::string>("roughness"),
        parameters,
        1,
        nullptr);

    return {length, diameter, casing_thickness, pipe_thickness, roughness};
}
}  // namespace WellboreSimulator
}  // namespace ProcessLib
