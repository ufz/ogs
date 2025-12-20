// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/ConfigTree.h"
#include "BishopsPowerLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<BishopsPowerLaw> createBishopsPowerLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "BishopsPowerLaw");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create BishopsPowerLaw property {:s}.", property_name);

    auto const exponent =
        //! \ogs_file_param{properties__property__BishopsPowerLaw__exponent}
        config.getConfigParameter<double>("exponent");

    return std::make_unique<MaterialPropertyLib::BishopsPowerLaw>(
        std::move(property_name), exponent);
}
}  // namespace MaterialPropertyLib
