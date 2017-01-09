/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

namespace MathLib
{

class LisMatrix;
class LisVector;

/**
 * \brief Linear solver using Lis (http://www.ssisc.org/lis/)
 *
 */
class LisLinearSolver final
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
    LisLinearSolver(const std::string solver_name = "",
                    BaseLib::ConfigTree const*const option = nullptr);

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LisOption &option) { _lis_option = option; }

    bool solve(LisMatrix& A, LisVector &b, LisVector &x);

private:
    LisOption _lis_option;
};

} // MathLib

#endif //LISLINEARSOLVER_H_

