/*!
   \file  PETScTools.cpp
   \brief Definitions of funstions related to PETSc solver interface.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013


  \copyright
   Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "PETScTools.h"

#include "logog/include/logog.hpp"

#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{

void  applyKnownSolution(PETScMatrix &A, PETScVector &b,  PETScVector &x,
                         const std::vector<int> &_vec_knownX_id,
                         const std::vector<double> &_vec_knownX_x)
{
   const PetscInt ni = static_cast<PetscInt> (_vec_knownX_id.size());


   A.zeroRows_in_Matrix(ni, &_vec_knownX_id[0]);
   A.finalAssemble();


   x.setValues(ni, &_vec_knownX_id[0], &_vec_knownX_x[0], INSERT_VALUES);
   b.setValues(ni, &_vec_knownX_id[0], &_vec_knownX_x[0], INSERT_VALUES);


   x.finalAssemble();
   b.finalAssemble();

}

} // end of namespace MathLib



