/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISLINEARSOLVER_H_
#define LISLINEARSOLVER_H_

#include <vector>
#include <string>

#include <lis.h>

#include "BaseLib/ConfigTree.h"

#include "LisOption.h"
#include "LisVector.h"
#include "LisMatrix.h"

namespace MathLib
{

/**
 * \brief Linear solver using Lis (http://www.ssisc.org/lis/)
 *
 */
class LisLinearSolver final
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
    LisLinearSolver(LisMatrix &A, const std::string solver_name = "",
                    BaseLib::ConfigTree const*const option = nullptr);

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LisOption &option) { _lis_option = option; }

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(LisVector &b, LisVector &x);


private:
    LisMatrix& _A;
    LisOption _lis_option;
};

} // MathLib

#endif //LISLINEARSOLVER_H_

