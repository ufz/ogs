/**
 * \file
 * \author Thomas Fischer
 * \date   2011-01-07
 * \brief  Definition of the DenseDirectLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef DENSEDIRECTLINEARSOLVER_H_
#define DENSEDIRECTLINEARSOLVER_H_

#include "DirectLinearSolver.h"

#include "../Dense/DenseMatrix.h"
#include "../Dense/DenseVector.h"

namespace MathLib {

/**
 * DenseDirectLinearSolver class provide interface to solve linear equations
 * represented with a dense matrix and vector. The class currently supports
 * only Gauss elimination algorithm.
 */
class DenseDirectLinearSolver: public MathLib::DirectLinearSolver {
public:
	DenseDirectLinearSolver() {};
	virtual ~DenseDirectLinearSolver() {};

    /**
     * solve a given system of linear equations
     *
     * @param A     Coefficient matrix
     * @param b     RHS vector
     * @param x     Solution vector
     */
    virtual void solve(DenseMatrix<double> &A, DenseVector<double> &b, DenseVector<double> &x);
};

}

#endif /* DENSEDIRECTLINEARSOLVER_H_ */
