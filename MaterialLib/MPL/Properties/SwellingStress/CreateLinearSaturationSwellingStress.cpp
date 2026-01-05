// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateLinearSaturationSwellingStress.h"

#include "BaseLib/ConfigTree.h"
#include "LinearSaturationSwellingStress.h"
#include "MaterialLib/MPL/Property.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createLinearSaturationSwellingStress(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "LinearSaturationSwellingStress");

    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create LinearSaturationSwellingStress phase property {:s}.",
         property_name);

    auto const coefficient =
        //! \ogs_file_param{properties__property__LinearSaturationSwellingStress__coefficient}
        config.getConfigParameter<double>("coefficient");

    auto const reference_saturation =
        //! \ogs_file_param{properties__property__LinearSaturationSwellingStress__reference_saturation}
        config.getConfigParameter<double>("reference_saturation");

    return std::make_unique<LinearSaturationSwellingStress>(
        property_name, coefficient, reference_saturation);
}
}  // namespace MaterialPropertyLib
