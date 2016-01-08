/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EIGENLINEARSOLVER_H_
#define EIGENLINEARSOLVER_H_

#include <vector>


#include "BaseLib/ConfigTreeNew.h"
#include "EigenVector.h"
#include "EigenOption.h"

namespace MathLib
{

class EigenMatrix;

class EigenLinearSolver final
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
    EigenLinearSolver(EigenMatrix &A, const std::string solver_name = "",
                      BaseLib::ConfigTreeNew const*const option = nullptr);

    ~EigenLinearSolver()
    {
        delete _solver;
    }

    /**
     * parse linear solvers configuration
     */
    void setOption(const BaseLib::ConfigTreeNew& option);

    /**
     * copy linear solvers options
     */
    void setOption(const EigenOption &option) { _option = option; }

    /**
     * get linear solver options
     */
    EigenOption &getOption() { return _option; }

    /**
     * solve a given linear equations
     *
     * @param b     RHS vector
     * @param x     Solution vector
     */
    void solve(EigenVector &b, EigenVector &x);

protected:
    class IEigenSolver
    {
    public:
        virtual ~IEigenSolver() = default;
        /**
         * execute a linear solver
         */
        virtual void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &) = 0;
    };

    EigenOption _option;
    IEigenSolver* _solver;
};

} // MathLib

#endif //EIGENLINEARSOLVER_H_

