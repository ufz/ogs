/**
 * \file
 * \author Norbert Grunwald
 * \date   Jan, 2021
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BaseLib/ConfigTree.h"
#include "IdealGasLawBinaryMixture.h"

namespace MaterialPropertyLib
{
std::unique_ptr<IdealGasLawBinaryMixture> createIdealGasLawBinaryMixture(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "IdealGasLawBinaryMixture");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create IdealGasLawBinaryMixture medium property {:s}.",
         property_name);

    //! \ogs_file_param_special{properties__property__IdealGasLawBinaryMixture}
    return std::make_unique<IdealGasLawBinaryMixture>(std::move(property_name));
}
}  // namespace MaterialPropertyLib
