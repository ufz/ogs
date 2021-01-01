/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "SaturationLiakopoulos");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create SaturationLiakopoulos medium property {:s}.", property_name);

    //! \ogs_file_param_special{properties__property__SaturationLiakopoulos}
    return std::make_unique<SaturationLiakopoulos>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
