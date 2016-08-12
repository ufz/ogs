/*!
   \file  SolidDensityType.h
   \brief Declaration of solid density types.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef SOLID_DENSITY_TYPE_H_
#define SOLID_DENSITY_TYPE_H_

namespace MaterialLib
{
enum class SolidDensityType
{
    CONSTANT = 0,
    LINEAR
};
}  // end namespace
#endif
