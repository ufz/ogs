// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
