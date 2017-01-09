/*!
   \file  PETScTools.cpp
   \brief Definition of a function related to PETSc solver interface to assign
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

#include "PETScTools.h"

namespace MathLib
{

void applyKnownSolution(PETScMatrix &A, PETScVector &b, PETScVector &x,
                        const std::vector<PetscInt>  &vec_knownX_id,
                        const std::vector<PetscScalar> &vec_knownX_x)
{
    A.finalizeAssembly();

    A.setRowsColumnsZero(vec_knownX_id);
    A.finalizeAssembly();

    x.finalizeAssembly();
    b.finalizeAssembly();
    if(vec_knownX_id.size() > 0)
    {
        x.set(vec_knownX_id, vec_knownX_x);
        b.set(vec_knownX_id, vec_knownX_x);
    }

    x.finalizeAssembly();
    b.finalizeAssembly();
}

} // end of namespace MathLib


