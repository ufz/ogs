/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef APPLICATIONS_NUMERICSCONFIG_H_
#define APPLICATIONS_NUMERICSCONFIG_H_

#include <type_traits>

/**
 * This file provides a configuration of the global matrix/vector and
 * corresponding linear solver, and the global executer types.
 * The configuration is collected in the GlobalSetupType being a particular
 * instantiation of the AssemblerLib::GlobalSetup template.
 * The existence of the GlobalSetupType is checked at the end of the file.
 */

//
// Global vector/matrix types and linear solver.
//
#ifdef USE_LIS

    #include "MathLib/LinAlg/Lis/LisMatrix.h"
    #include "MathLib/LinAlg/Lis/LisVector.h"
    #include "MathLib/LinAlg/Lis/LisTools.h"
    #include "MathLib/LinAlg/Lis/LisLinearSolver.h"

namespace detail
{
    using GlobalVectorType = MathLib::LisVector;
    using GlobalMatrixType = MathLib::LisMatrix;
    using GlobalMatrixAndNodeAdjacencyTableBuilder
          = MathLib::LisMatrixAndNodeAdjacencyTableBuilder;
    using LinearSolverType = MathLib::LisLinearSolver;
}

#else    // USE_LIS
#ifdef OGS_USE_EIGEN
    #include "MathLib/LinAlg/Eigen/EigenVector.h"
    #include "MathLib/LinAlg/Eigen/EigenMatrix.h"
    #include "MathLib/LinAlg/Eigen/EigenTools.h"
    #include "MathLib/LinAlg/Eigen/EigenLinearSolver.h"
namespace detail
{
    using GlobalVectorType = MathLib::EigenVector;
    using GlobalMatrixType = MathLib::EigenMatrix;
    using GlobalMatrixAndNodeAdjacencyTableBuilder
          = MathLib::EigenMatrixAndNodeAdjacencyTableBuilder;
    using LinearSolverType = MathLib::EigenLinearSolver;
}
#else   // OGS_USE_EIGEN
    #include "MathLib/LinAlg/Dense/DenseVector.h"
    #include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
    #include "MathLib/LinAlg/Dense/DenseTools.h"
    #include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
namespace detail
{
    using GlobalVectorType = MathLib::DenseVector<double>;
    using GlobalMatrixType = MathLib::GlobalDenseMatrix<double>;
    using GlobalMatrixAndNodeAdjacencyTableBuilder
          = MathLib::DenseMatrixAndNodeAdjacencyTableBuilder;
    using LinearSolverType =
        MathLib::GaussAlgorithm<GlobalMatrixType, GlobalVectorType>;
}

#endif    // USE_LIS
#endif    // OGS_USE_EIGEN


//
// Global vector/matrix builder.
//

#include "AssemblerLib/VectorMatrixBuilder.h"
namespace detail
{
using GlobalVectorMatrixBuilderType =
        AssemblerLib::VectorMatrixBuilder<
            GlobalMatrixType,
            GlobalMatrixAndNodeAdjacencyTableBuilder,
            GlobalVectorType>;
}

//
// Global executor
//
#include "AssemblerLib/SerialExecutor.h"
namespace detail
{
using GlobalExecutorType = AssemblerLib::SerialExecutor;
}

///
/// Global setup collects the previous configuration in single place.
///
#include "AssemblerLib/GlobalSetup.h"
using GlobalSetupType =
    AssemblerLib::GlobalSetup<
        detail::GlobalVectorMatrixBuilderType,
        detail::GlobalExecutorType,
        detail::LinearSolverType>;


//
// Check the configuration
//
static_assert(std::is_class<GlobalSetupType>::value,
              "GlobalSetupType was not defined.");
static_assert(std::is_integral<detail::GlobalMatrixType::IndexType>::value,
              "The index type for global matrices is not an integral type.");
static_assert(std::is_integral<detail::GlobalVectorType::IndexType>::value,
              "The index type for global vectors is not an integral type.");
static_assert(std::is_same<detail::GlobalMatrixType::IndexType,
                           detail::GlobalVectorType::IndexType>::value,
              "The global matrix and vector index types do not match.");
// Both types are integral types and equal, define a single GlobalIndexType.

/// A type used for indexing of global vectors and matrices. It is equal to the
/// GlobalMatrixType::IndexType and the GlobalVectorType::IndexType.
using GlobalIndexType = detail::GlobalMatrixType::IndexType;

#endif  // APPLICATIONS_NUMERICSCONFIG_H_
