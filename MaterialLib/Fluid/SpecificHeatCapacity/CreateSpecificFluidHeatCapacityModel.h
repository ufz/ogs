/**
 *  \brief A function for creating a specific heat capacity model for fluid
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file CreateSpecificFluidHeatCapacityModel.h
 *
 */

#ifndef OGS_CREATESPECIFICFLUIDHEATCAPACITYMODEL_H
#define OGS_CREATESPECIFICFLUIDHEATCAPACITYMODEL_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace Fluid
{
class FluidProperty;

/**
 *  Create a specific heat capacity model
 *  \param config  ConfigTree object has a tag of `<specific_heat_capacity>`
 */
std::unique_ptr<FluidProperty> createSpecificFluidHeatCapacityModel(
    BaseLib::ConfigTree const& config);

}  // end namespace
}  // end namespace

#endif /* OGS_CREATESPECIFICFLUIDHEATCAPACITY_H */
