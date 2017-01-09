/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "SparsityPattern.h"

//
// Global vector/matrix types and linear solver.
//
#if defined(USE_LIS)

    #include "MathLib/LinAlg/Eigen/EigenMatrix.h"
    #include "MathLib/LinAlg/Eigen/EigenVector.h"
    #include "MathLib/LinAlg/EigenLis/EigenLisLinearSolver.h"

    using GlobalVector = MathLib::EigenVector;
    using GlobalMatrix = MathLib::EigenMatrix;

    using GlobalLinearSolver = MathLib::EigenLisLinearSolver;

#elif defined(USE_PETSC)
    #include "MathLib/LinAlg/PETSc/PETScVector.h"
    #include "MathLib/LinAlg/PETSc/PETScMatrix.h"
    #include "MathLib/LinAlg/PETSc/PETScLinearSolver.h"

    using GlobalVector = MathLib::PETScVector;
    using GlobalMatrix = MathLib::PETScMatrix;

    using GlobalLinearSolver = MathLib::PETScLinearSolver;

#elif defined(OGS_USE_EIGEN)
    #include "MathLib/LinAlg/Eigen/EigenVector.h"
    #include "MathLib/LinAlg/Eigen/EigenMatrix.h"
    #include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"

    using GlobalVector = MathLib::EigenVector;
    using GlobalMatrix = MathLib::EigenMatrix;

    using GlobalLinearSolver = MathLib::EigenLinearSolver;

#endif    // OGS_USE_EIGEN


/// A type used for indexing of global vectors and matrices. It is equal to the
/// GlobalMatrix::IndexType and the GlobalVector::IndexType.
static_assert(std::is_integral<GlobalMatrix::IndexType>::value,
              "The index type for global matrices is not an integral type.");
static_assert(std::is_integral<GlobalVector::IndexType>::value,
              "The index type for global vectors is not an integral type.");
static_assert(std::is_same<GlobalMatrix::IndexType,
                           GlobalVector::IndexType>::value,
              "The global matrix and vector index types do not match.");
// Both types are integral types and equal, define a single GlobalIndexType.
using GlobalIndexType = GlobalMatrix::IndexType;

using GlobalSparsityPattern = MathLib::SparsityPattern<GlobalIndexType>;
