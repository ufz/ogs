/*!
   \file  PermeabilityType.h
   \brief Declaration of permeability types.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef PERMEABILITY_TYPE_H_
#define PERMEABILITY_TYPE_H_

namespace MaterialLib
{

/// Intrinsic permeability or hydraulic conductivity type
enum class PermeabilityType
{
    ISOTROPIC = 0,
    ANISOTROPIC = 1
};

} // end namespace
#endif
