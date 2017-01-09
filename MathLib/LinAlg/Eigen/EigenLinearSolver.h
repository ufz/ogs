/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENLINEARSOLVER_H_
#define EIGENLINEARSOLVER_H_

#include <vector>

#include "BaseLib/ConfigTree.h"
#include "EigenOption.h"

namespace MathLib
{

class EigenMatrix;
class EigenVector;

class EigenLinearSolverBase;

class EigenLinearSolver final
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
    EigenLinearSolver(const std::string& solver_name,
                      BaseLib::ConfigTree const*const option);

    ~EigenLinearSolver();

    /**
     * parse linear solvers configuration
     */
    void setOption(const BaseLib::ConfigTree& option);

    /**
     * copy linear solvers options
     */
    void setOption(const EigenOption &option) { _option = option; }

    /**
     * get linear solver options
     */
    EigenOption &getOption() { return _option; }

    bool solve(EigenMatrix &A, EigenVector& b, EigenVector &x);

protected:
    EigenOption _option;
    std::unique_ptr<EigenLinearSolverBase> _solver;
};

} // MathLib

#endif //EIGENLINEARSOLVER_H_

