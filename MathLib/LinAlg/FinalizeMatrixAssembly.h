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

#ifndef FINALIZEMATRIXASSEMBLY_H_
#define FINALIZEMATRIXASSEMBLY_H_

namespace MathLib
{

template <typename MAT_T>
bool finalizeMatrixAssembly(MAT_T &)
{
	return true;
}

} // MathLib

#endif /* FINALIZEMATRIXASSEMBLY_H_ */
