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
#include "RelPermLiakopoulos.h"

namespace MaterialPropertyLib
{
std::unique_ptr<RelPermLiakopoulos> createRelPermLiakopoulos(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "RelPermLiakopoulos");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create RelPermLiakopoulos medium property {:s}.", property_name);

    //! \ogs_file_param_special{properties__property__RelPermLiakopoulos}
    return std::make_unique<RelPermLiakopoulos>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
