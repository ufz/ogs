/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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
    #include "MathLib/LinAlg/Lis/LisLinearSolver.h"

namespace detail
{
    using GlobalVectorType = MathLib::LisVector;
    using GlobalMatrixType = MathLib::LisMatrix;

    using LinearSolverType = MathLib::LisLinearSolver;
}

#else    // USE_LIS
    #include "MathLib/LinAlg/Dense/DenseVector.h"
    #include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
    #include "MathLib/LinAlg/Solvers/GaussAlgorithm.h"
namespace detail
{
    using GlobalVectorType = MathLib::DenseVector<double>;
    using GlobalMatrixType = MathLib::GlobalDenseMatrix<double>;

    using LinearSolverType =
        MathLib::GaussAlgorithm<GlobalMatrixType, GlobalVectorType>;
}

#endif    // USE_LIS


//
// Global executor and vector/matrix builder.
//
#include "AssemblerLib/SerialExecutor.h"
#include "AssemblerLib/SerialVectorMatrixBuilder.h"

namespace detail
{
using GlobalExecutorType = AssemblerLib::SerialExecutor;

using GlobalVectorMatrixBuilderType =
        AssemblerLib::SerialVectorMatrixBuilder<
            GlobalMatrixType,
            GlobalVectorType>;
}


#include "AssemblerLib/GlobalSetup.h"

///
/// Global setup collects the previous configuration in single place.
///
using GlobalSetupType =
    AssemblerLib::GlobalSetup<
        detail::GlobalVectorMatrixBuilderType,
        detail::GlobalExecutorType>;


//
// Check the configuration
//
static_assert(std::is_class<GlobalSetupType>::value,
        "GlobalSetupType was not defined.");
#endif  // APPLICATIONS_NUMERICSCONFIG_H_
