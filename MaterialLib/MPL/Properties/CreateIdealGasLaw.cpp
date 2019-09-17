/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "IdealGasLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<IdealGasLaw> createIdealGasLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "IdealGasLaw");
    DBUG("Create IdealGasLaw medium property");
    return std::make_unique<IdealGasLaw>();
}
}  // namespace MaterialPropertyLib
