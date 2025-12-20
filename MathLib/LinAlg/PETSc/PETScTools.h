// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "MathLib/LinAlg/LinAlgEnums.h"
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
   \param mode           Provided for compatibility with Eigen, unused
*/
void applyKnownSolution(PETScMatrix& A, PETScVector& b, PETScVector& x,
                        const std::vector<PetscInt>& vec_knownX_id,
                        const std::vector<PetscScalar>& vec_knownX_x,
                        DirichletBCApplicationMode const mode);
}  // end of namespace MathLib
