/**
 * \file
 * \author Norbert Grunwald
 * \date   Sep 10, 2019
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "Constant.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Constant> createConstant(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Constant");
    DBUG("Create Constant property");
    std::vector<double> const values =
        //! \ogs_file_param{properties__property__Constant__value}
        config.getConfigParameter<std::vector<double>>("value");

    return std::make_unique<Constant>(fromVector(values));
}
}  // namespace MaterialPropertyLib
