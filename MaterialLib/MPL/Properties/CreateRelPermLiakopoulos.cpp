/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
    // check is reading the parameter, not peeking it...
    //! \ogs_file_param{prj__media__medium__properties__property__RelPermLiakopoulos}
    config.checkConfigParameter("type", "RelPermLiakopoulos");
    DBUG("Create RelPermLiakopoulos medium property");

    return std::make_unique<RelPermLiakopoulos>();
}
}  // namespace MaterialPropertyLib
