/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "BishopsPowerLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<BishopsPowerLaw> createBishopsPowerLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "BishopsPowerLaw");
    DBUG("Create BishopsPowerLaw property");

    auto const exponent =
        //! \ogs_file_param{properties__property__BishopsPowerLaw__exponent}
        config.getConfigParameter<double>("exponent");

    return std::make_unique<MaterialPropertyLib::BishopsPowerLaw>(exponent);
}
}  // namespace MaterialPropertyLib
