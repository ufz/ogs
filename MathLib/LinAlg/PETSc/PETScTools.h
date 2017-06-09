/*!
   \file  PETScTools.h
   \brief Declaration of a function related to PETSc solver interface to assign
         the Dirichlet boundary conditions.

   \author Wenqing Wang
   \version
   \date Nov 2011 - Sep 2013

   \copyright
    Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#pragma once

#include <vector>

#include "PETScMatrix.h"
#include "PETScVector.h"

namespace MathLib
{
/*!
   \brief apply known solutions to a system of linear equations

   \param A              Coefficient matrix
   \param b              RHS vector
   \param x              Solution vector
   \param vec_knownX_id  A vector of known solution entry IDs
   \param vec_knownX_x   A vector of known solutions
*/
void applyKnownSolution(PETScMatrix& A, PETScVector& b, PETScVector& x,
                        const std::vector<PetscInt>& vec_knownX_id,
                        const std::vector<PetscScalar>& vec_knownX_x);
}  // end of namespace MathLib
