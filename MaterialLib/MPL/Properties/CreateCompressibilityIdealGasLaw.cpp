/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *  Created on January 28, 2020, 16:05 PM
 */

#include "BaseLib/ConfigTree.h"
#include "CompressibilityIdealGasLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<CompressibilityIdealGasLaw> createCompressibilityIdealGasLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "CompressibilityIdealGasLaw");
    DBUG("Create CompressibilityIdealGasLaw medium property");
    return std::make_unique<CompressibilityIdealGasLaw>();
}
}  // namespace MaterialPropertyLib
