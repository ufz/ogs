/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BaseLib/ConfigTree.h"
#include "SaturationLiakopoulos.h"

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationLiakopoulos> createSaturationLiakopoulos(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__media__medium__properties__property__SaturationLiakopoulos}
    config.checkConfigParameter("type", "SaturationLiakopoulos");
    DBUG("Create SaturationLiakopoulos medium property");

    return std::make_unique<SaturationLiakopoulos>();
}
}  // namespace MaterialPropertyLib
