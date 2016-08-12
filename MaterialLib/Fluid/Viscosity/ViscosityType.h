/*!
   \file  ViscosityType.h
   \brief Declaration of viscosity types.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef VISCOSITY_TYPE_H_
#define VISCOSITY_TYPE_H_

namespace MaterialLib
{
namespace Fluid
{
enum class ViscosityType
{
    CONSTANT = 0,
    LINEAR_PRESSURE_DEPENDENT,
    TEMPERATURE_DEPENDENT,
    VOGEL
};

}  // end namespace
}  // end namespace
#endif
