/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateConstantPermeability.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "ConstantPermeability.h"

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createConstantPermeability(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fracture_properties__permeability_model__type}
    config.checkConfigParameter("type", "ConstantPermeability");

    auto const permeability =
        //! \ogs_file_param{material__fracture_properties__permeability_model__ConstantPermeability__value}
        config.getConfigParameter<double>("value");

    if (permeability < 0)
    {
        OGS_FATAL(
            "The permeability parameter must be non-negative. Given value "
            "{:g}.",
            permeability);
    }

    return std::make_unique<ConstantPermeability>(permeability);
}
}  // namespace MaterialLib::Fracture::Permeability
