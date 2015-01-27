/*!
   \file  DensityType.h
   \brief Declaration of density types.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef DENSITY_TYPE_H_
#define DENSITY_TYPE_H_

namespace MaterialLib
{

enum class DensityType
{
    CONSTANT,
    SOLID_LINEAR,
    IDEAL_GAS
};

} // end namespace
#endif
