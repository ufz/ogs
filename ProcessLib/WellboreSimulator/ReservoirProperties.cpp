// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ReservoirProperties.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
namespace WellboreSimulator
{
ReservoirProperties createReservoirProperties(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    auto const& temperature = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties__temperature}
        config.getConfigParameter<std::string>("temperature"),
        parameters,
        1,
        nullptr);

    auto const& thermal_conductivity = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties__thermal_conductivity}
        config.getConfigParameter<std::string>("thermal_conductivity"),
        parameters,
        1,
        nullptr);

    auto const& density = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties__density}
        config.getConfigParameter<std::string>("density"),
        parameters,
        1,
        nullptr);

    auto const& specific_heat_capacity = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties__specific_heat_capacity}
        config.getConfigParameter<std::string>("specific_heat_capacity"),
        parameters,
        1,
        nullptr);

    auto const& pressure = ParameterLib::findParameter<double>(
        //! \ogs_file_param{prj__processes__process__WELLBORE_SIMULATOR__reservoir_properties__pressure}
        config.getConfigParameter<std::string>("pressure"),
        parameters,
        1,
        nullptr);

    return {temperature, thermal_conductivity, density, specific_heat_capacity,
            pressure};
}
}  // namespace WellboreSimulator
}  // namespace ProcessLib
