/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief  Definition of the LisLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef LISLINEARSOLVER_H_
#define LISLINEARSOLVER_H_

#include <vector>
#include <boost/property_tree/ptree.hpp>

#include "lis.h"

#include "LisOption.h"
#include "LisVector.h"
#include "LisMatrix.h"

namespace MathLib
{

/**
 * \brief Linear solver using Lis (http://www.ssisc.org/lis/)
 *
 */
class LisLinearSolver
{
public:
    LisLinearSolver();

    /**
     *
     */
    virtual ~LisLinearSolver();

    /**
     * configure linear solvers
     * @param option
     */
    virtual void setOption(const boost::property_tree::ptree &option);

    /// solve a given linear equations
    virtual void solve(LisMatrix &A, LisVector &b, LisVector &x);

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LisOption &option) { _option = option; }

    /**
     * get linear solver options
     * @return
     */
    LisOption &getOption() { return _option; }

private:
    LisOption _option;
};

} // MathLib

#endif //LISLINEARSOLVER_H_

