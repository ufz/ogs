/*!
   \file  FluidDensityType.h
   \brief Declaration of fluid density types.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef FLUIDDENSITY_TYPE_H_
#define FLUIDDENSITY_TYPE_H_

namespace MaterialLib
{
namespace Fluid
{
enum class FluidDensityType
{
    CONSTANT = 0,
    IDEAL_GAS,
    LIQUID_DENSITY,
    LINEAR_TEMPERATURE_DEPENDENT
};
}
}  // end namespace
#endif
