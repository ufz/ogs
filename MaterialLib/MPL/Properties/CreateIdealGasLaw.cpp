// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "BaseLib/ConfigTree.h"
#include "IdealGasLaw.h"

namespace MaterialPropertyLib
{
std::unique_ptr<IdealGasLaw> createIdealGasLaw(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "IdealGasLaw");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create IdealGasLaw medium property {:s}.", property_name);
    //! \ogs_file_param_special{properties__property__IdealGasLaw}
    return std::make_unique<IdealGasLaw>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
