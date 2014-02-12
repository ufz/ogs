/*!

  \file FinalizeVectorAssembly.h
  \author Wenqing Wang
  \date Oct, 2013

   \copyright
    Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/LICENSE.txt
*/

#ifndef FINALIZEVECTORASSEMBLY_H_
#define FINALIZEVECTORASSEMBLY_H_

namespace MathLib
{

template <typename VEC_T>
void  finalizeVectorAssembly(VEC_T & vec)
{
    vec.finalizeAssembly();
}

} // MathLib

#endif
