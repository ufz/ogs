/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateCubicLaw.h"

#include "BaseLib/ConfigTree.h"
#include "CubicLaw.h"

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createCubicLaw(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__fracture_properties__permeability_model__type}
    config.checkConfigParameter("type", "CubicLaw");

    //! \ogs_file_param_special{material__fracture_properties__permeability_model__CubicLaw}
    return std::make_unique<CubicLaw>();
}
}  // namespace MaterialLib::Fracture::Permeability
