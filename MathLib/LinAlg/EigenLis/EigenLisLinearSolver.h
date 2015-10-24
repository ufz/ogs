/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENLISLINEARSOLVER_H_
#define EIGENLISLINEARSOLVER_H_

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
     * @param A           Coefficient matrix object
     * @param solver_name A name used as a prefix for command line options
     *                    if there are such options available.
     * @param option      A pointer to a linear solver option. In case you omit
     *                    this argument, default settings follow those of
     *                    LisOption struct.
     */
    EigenLisLinearSolver(EigenMatrix &A, const std::string solver_name = "",
                         BaseLib::ConfigTree const*const option = nullptr);

    /**
     * parse linear solvers configuration
     */
    void setOption(BaseLib::ConfigTree const& option);

    /**
     * copy linear solvers options
     */
    void setOption(const LisOption &option) { _option = option; }

    /**
     * get linear solver options
     */
    LisOption &getOption() { return _option; }

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(EigenVector &b, EigenVector &x);

private:
    EigenMatrix& _A;
    LisOption _option;
};

} // MathLib

#endif //EIGENLISLINEARSOLVER_H_

