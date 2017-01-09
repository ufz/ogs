/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_NUMERICSCONFIG_H_
#define NUMLIB_NUMERICSCONFIG_H_

/**
 * This file provides a configuration of the global matrix/vector and
 * corresponding linear solver, and the global executer types.
 */

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

//
// Global executor
//
#include "NumLib/Assembler/SerialExecutor.h"
using GlobalExecutor = NumLib::SerialExecutor;



#endif  // NUMLIB_NUMERICSCONFIG_H_
