/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef APPLICATIONS_NUMERICSCONFIG_H_
#define APPLICATIONS_NUMERICSCONFIG_H_

#include <type_traits>
#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

/**
 * This file provides a configuration of the global matrix/vector and
 * corresponding linear solver, and the global executer types.
 * The configuration is collected in the GlobalSetupType being a particular
 * instantiation of the ProcessLib::GlobalSetup template.
 * The existence of the GlobalSetupType is checked at the end of the file.
 */

//
// Global vector/matrix builder.
//

#include "NumLib/Assembler/VectorMatrixBuilder.h"
namespace detail
{
using GlobalVectorMatrixBuilderType =
        NumLib::VectorMatrixBuilder<
            GlobalMatrixType,
            GlobalVectorType>;
}

//
// Global executor
//
#include "NumLib/Assembler/SerialExecutor.h"
namespace detail
{
using GlobalExecutorType = NumLib::SerialExecutor;
}

///
/// Global setup collects the previous configuration in single place.
///
#include "GlobalSetup.h"
using GlobalSetupType =
    NumLib::GlobalSetup<
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

#endif  // APPLICATIONS_NUMERICSCONFIG_H_
