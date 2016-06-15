/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NUMERICSCONFIG_H_
#define NUMLIB_NUMERICSCONFIG_H_

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

/**
 * This file provides a configuration of the global matrix/vector and
 * corresponding linear solver, and the global executer types.
 * The configuration is collected in the GlobalSetupType being a particular
 * instantiation of the ProcessLib::GlobalSetup template.
 * The existence of the GlobalSetupType is checked at the end of the file.
 */

//
// Global executor
//
#include "NumLib/Assembler/SerialExecutor.h"
namespace detail
{
using GlobalExecutorType = NumLib::SerialExecutor;
}

#endif  // NUMLIB_NUMERICSCONFIG_H_
