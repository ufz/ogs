/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include <lis.h>

#include "BaseLib/ConfigTree.h"
#include "MathLib/LinAlg/Lis/LisOption.h"

namespace MathLib
{

class EigenVector;
class EigenMatrix;

/**
 * Linear solver using Lis library with Eigen matrix and vector objects
 */
class EigenLisLinearSolver final
{
public:
    /**
     * Constructor
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option      A pointer to a linear solver option. In case you omit
     *                    this argument, default settings follow those of
     *                    LisOption struct.
     */
    EigenLisLinearSolver(const std::string solver_name,
                         BaseLib::ConfigTree const*const option);

    /**
     * copy linear solvers options
     */
    void setOption(const LisOption &option) { _lis_option = option; }

    bool solve(EigenMatrix &A, EigenVector& b, EigenVector &x);

private:
    LisOption _lis_option;
};

} // MathLib
