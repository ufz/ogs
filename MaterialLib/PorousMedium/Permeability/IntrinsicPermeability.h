/*!
   \file  IntrisincPermeability.h
   \brief Declaration of class IntrisincPermeability.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#ifndef INTRINSIC_PERMEABILITY_H_
#define INTRINSIC_PERMEABILITY_H_

#include "MaterialLib/ConstantTensor.h"

namespace MaterialLib
{
namespace PorousMedia
{
template <typename T_MATRIX>
using IntrinsicPermeability = ConstantTensor<T_MATRIX>;

}  // end namespace
}  // end namespace
#endif
