/*!
   \file  IntrisincPermeability.h
   \brief Declaration of class IntrisincPermeability.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef INTRISINC_PERMEABILITY_H_
#define INTRISINC_PERMEABILITY_H_

#include <type_traits>

#include "MaterialLib/ConstantTensor.h"

namespace MaterialLib
{
template<typename T_MATRIX > using IntrisincPermeability
                                    =  ConstantTensor<T_MATRIX>;

} // end namespace
#endif
