/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 27, 2020, 12:27 PM
 */

#include "CreateCapillaryPressureLiakopoulos.h"

#include "BaseLib/ConfigTree.h"
#include "CapillaryPressureLiakopoulos.h"

namespace MaterialPropertyLib
{
std::unique_ptr<CapillaryPressureLiakopoulos>
createCapillaryPressureLiakopoulos(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CapillaryPressureLiakopoulos");
    DBUG("Create CapillaryPressureLiakopoulos medium property");

    //! \ogs_file_param_special{properties__property__CapillaryPressureLiakopoulos}
    return std::make_unique<CapillaryPressureLiakopoulos>();
}
}  // namespace MaterialPropertyLib
