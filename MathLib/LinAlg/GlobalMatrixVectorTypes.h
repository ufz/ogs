/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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

using GlobalVector = MathLib::EigenVector;
using GlobalMatrix = MathLib::EigenMatrix;

#elif defined(USE_PETSC)
#include "MathLib/LinAlg/PETSc/PETScMatrix.h"
#include "MathLib/LinAlg/PETSc/PETScVector.h"

using GlobalVector = MathLib::PETScVector;
using GlobalMatrix = MathLib::PETScMatrix;

#else
#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"

using GlobalVector = MathLib::EigenVector;
using GlobalMatrix = MathLib::EigenMatrix;

#endif

/// A type used for indexing of global vectors and matrices. It is equal to the
/// GlobalMatrix::IndexType and the GlobalVector::IndexType.
static_assert(std::is_integral_v<GlobalMatrix::IndexType>,
              "The index type for global matrices is not an integral type.");
static_assert(std::is_integral_v<GlobalVector::IndexType>,
              "The index type for global vectors is not an integral type.");
static_assert(std::is_same_v<GlobalMatrix::IndexType, GlobalVector::IndexType>,
              "The global matrix and vector index types do not match.");
// Both types are integral types and equal, define a single GlobalIndexType.
using GlobalIndexType = GlobalMatrix::IndexType;

using GlobalSparsityPattern = MathLib::SparsityPattern<GlobalIndexType>;
