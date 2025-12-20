// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
