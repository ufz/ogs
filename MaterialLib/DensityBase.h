/*!
   \file  DensityBase.h
   \brief Declaration of root class for material density.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef DENSITY_BASE_H_
#define DENSITY_BASE_H_

#include<memory>
#include<string>

namespace MaterialLib
{

enum class DensityType
{
    CONSTANT,
    SOLID_LINEAR,
    IDEAL_GAS
};

/// Base density class.
class DensityBase
{
      public:
        double getDensity() {return 0.;}
};

} // end namespace
#endif
