/**
 * \file
 * \author Norbert Grunwald
 * \date   07.09.2017
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Medium.h"

#include <string>
#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Parameter.h"

#include "Properties/Properties.h"

#include "CreatePhase.h"
#include "CreateProperty.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Medium> createMedium(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    // Parsing the phases
    // Properties of phases may be not required in all the cases.
    auto&& phases =
        //! \ogs_file_param{prj__media__medium__phases}
        createPhases(config.getConfigSubtreeOptional("phases"), parameters);

    // Parsing medium properties, overwriting the defaults.
    auto&& properties =
        //! \ogs_file_param{prj__media__medium__properties}
        createProperties(config.getConfigSubtreeOptional("properties"),
                         parameters);

    if (phases.empty() && !properties)
    {
        OGS_FATAL("Neither tag <phases> nor tag <properties> has been found.");
    }

    return std::make_unique<Medium>(std::move(phases), std::move(properties));
}

}  // namespace MaterialPropertyLib
