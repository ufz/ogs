/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 8:47 AM
 */

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

    return std::make_unique<LinearSaturationSwellingStress>(property_name,
                                                            coefficient);
}
}  // namespace MaterialPropertyLib
