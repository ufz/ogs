/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file LisLinearEquation.h
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#ifndef LISLINEAREQUATION_H_
#define LISLINEAREQUATION_H_

#include <string>

#include "lis.h"

#include "BaseLib/CodingTools.h"
#include "AbstractCRSLinearEquation.h"
#include "LIS_option.h"

namespace MathLib
{

/**
 * \brief Linear equation using Lis solver (http://www.ssisc.org/lis/) with CRS matrix
 *
 * This class does not support MPI version of Lis solver.
 */
class LisLinearEquation : public AbstractCRSLinearEquation<signed>
{
public:
    /**
     * initialize any internal library
     */
    static void initialize();

    /**
     * finalize any internal library
     */
    static void finalize();

    /**
     *
     */
    LisLinearEquation(size_t length, RowMajorSparsity* sp);

    /**
     *
     */
    virtual ~LisLinearEquation();

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const BaseLib::Options &option);

    /**
     * configure linear solvers
     * @param option
     */
    void setOption(const LIS_option &option)
    {
        _option = option;
    }

    /**
     * get linear solver options
     * @return
     */
    LIS_option &getOption() 
    {
        return _option;
    }

    //void gatherX(std::vector<double> &x);

protected:
    void solveEqs(CRSMatrix<double, signed> *A, double *rhs, double *x);

private:
    LIS_option _option;
    LIS_VECTOR bb,xx;
};


} // MathLib

#endif //LISLINEAREQUATION_H_

