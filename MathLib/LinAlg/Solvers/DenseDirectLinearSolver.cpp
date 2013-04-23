/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-18
 * \brief  Implementation of the DenseDirectLinearSolver class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DenseDirectLinearSolver.h"

#include "GaussAlgorithm.h"

namespace MathLib
{

void DenseDirectLinearSolver::solve(Matrix<double> &A, Vector<double> &b, Vector<double> &x)
{
    GaussAlgorithm gauss(A);
    x = b;
    gauss.execute(x.getData());
}

}

