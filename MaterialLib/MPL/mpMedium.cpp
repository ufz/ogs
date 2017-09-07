/**
 * \author Norbert Grunwald
 * \date   07.09.2017
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "mpMedium.h"

namespace MaterialPropertyLib
{

Medium::Medium(BaseLib::ConfigTree const& config)
{
    // A Medium consists of phases and properties only.
    // Parse the phase configurations and push them into the
    // Medium::_phases attribute;
    auto const phases_config = config.getConfigSubtree("phases");
    createPhases(phases_config);

    // Parse the property configurations. These are medium properties only.
    // The properties of the phases are handled in the respective phases
    auto const properties_config = config.getConfigSubtree("properties");
    createProperties(properties_config);
}

void Medium::createPhases(BaseLib::ConfigTree const& config)
{

}

void Medium::createProperties(BaseLib::ConfigTree const& config)
{

}

} // MaterialPropertyLib




