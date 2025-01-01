/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreatePengRobinson.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"
#include "PengRobinson.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPengRobinson(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "PengRobinson");
    DBUG("Create PengRobinson EOS.");

    auto const Tc =
        //! \ogs_file_param{properties__property__PengRobinson__critical_temperature}
        config.getConfigParameter<double>("critical_temperature");

    auto const pc =
        //! \ogs_file_param{properties__property__PengRobinson__critical_pressure}
        config.getConfigParameter<double>("critical_pressure");

    auto const omega =
        //! \ogs_file_param{properties__property__PengRobinson__acentric_factor}
        config.getConfigParameter<double>("acentric_factor");

    return std::make_unique<PengRobinson>(Tc, pc, omega);
}
}  // namespace MaterialPropertyLib
