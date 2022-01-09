/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreatePermeabilityModel.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "CreateConstantPermeability.h"
#include "CreateCubicLaw.h"

namespace MaterialLib::Fracture::Permeability
{
std::unique_ptr<Permeability> createPermeabilityModel(
    BaseLib::ConfigTree const& config)
{
    auto const permeability_model_type =
        //! \ogs_file_param{material__fracture_properties__permeability_model__type}
        config.peekConfigParameter<std::string>("type");
    if (permeability_model_type == "ConstantPermeability")
    {
        return MaterialLib::Fracture::Permeability::createConstantPermeability(
            config);
    }
    if (permeability_model_type == "CubicLaw")
    {
        return MaterialLib::Fracture::Permeability::createCubicLaw(config);
    }
    OGS_FATAL("Unknown fracture permeability model type \"{:s}\".",
              permeability_model_type);
}
}  // namespace MaterialLib::Fracture::Permeability
