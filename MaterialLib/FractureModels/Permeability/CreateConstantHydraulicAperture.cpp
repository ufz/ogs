/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateConstantHydraulicAperture.h"
#include "BaseLib/ConfigTree.h"

#include "ConstantHydraulicAperture.h"

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createConstantHydraulicAperture(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{prj__processes__process__type}
    config.checkConfigParameter("type", "ConstantHydraulicAperture");

    return std::make_unique<ConstantHydraulicAperture>();
}
}  // namespace MaterialLib::Fracture::Permeability
