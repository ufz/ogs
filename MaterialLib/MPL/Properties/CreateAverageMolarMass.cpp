/**
 * \file
 * \author Norbert Grunwald
 * \date   Jul 07 2020
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "AverageMolarMass.h"
#include "BaseLib/ConfigTree.h"

namespace MaterialPropertyLib
{
std::unique_ptr<AverageMolarMass> createAverageMolarMass(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "AverageMolarMass");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create AverageMolarMass medium property");
    //! \ogs_file_param_special{properties__property__AverageMolarMass}
    return std::make_unique<AverageMolarMass>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
