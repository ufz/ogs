/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "BishopsSaturationCutoff.h"

namespace MaterialPropertyLib
{
std::unique_ptr<BishopsSaturationCutoff> createBishopsSaturationCutoff(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "BishopsSaturationCutoff");
    DBUG("Create BishopsSaturationCutoff property");

    auto const cutoff_value =
        //! \ogs_file_param{properties__property__BishopsSaturationCutoff__cutoff_value}
        config.getConfigParameter<double>("cutoff_value");

    return std::make_unique<MaterialPropertyLib::BishopsSaturationCutoff>(
        cutoff_value);
}
}  // namespace MaterialPropertyLib
