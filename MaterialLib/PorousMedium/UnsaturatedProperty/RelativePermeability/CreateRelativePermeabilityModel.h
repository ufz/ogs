/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateRelativePermeabilityModel.h
 *
 * Created on November 2, 2016, 11:43 AM
 */

#ifndef OGS_CREATE_RELATIVE_PERMEABILITY_MODEL_H
#define OGS_CREATE_RELATIVE_PERMEABILITY_MODEL_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace PorousMedium
{
class RelativePermeability;
/// Create a capillary pressure model
/// \param config  ConfigTree object has a tag of <relative_permeability>
std::unique_ptr<RelativePermeability> createRelativePermeabilityModel(
    BaseLib::ConfigTree const& config);
}
}  // end namespace

#endif /* OGS_CREATE_RELATIVE_PERMEABILITY_MODEL_H */
