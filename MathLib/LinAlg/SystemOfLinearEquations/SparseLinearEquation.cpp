/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file SparseLinearEquation.cpp
 *
 * Created on 2012-06-25 by Norihiro Watanabe
 */

#include "SparseLinearEquation.h"

#include <iostream>
#include <algorithm>

#include "MathLib/LinAlg/Solvers/CG.h"
#include "MathLib/LinAlg/Solvers/BiCGStab.h"
#include "MathLib/LinAlg/Sparse/sparse_io.h"

namespace MathLib
{


void SparseLinearEquation::setOption(const BaseLib::Options &option)
{
    const BaseLib::Options *op = option.getSubGroup("LinearSolver");
    if (op==0) {
        return;
    }

    if (op->hasOption("solver_type"))
        _option.solver_type = _option.getSolverType(op->getOption("solver_type"));
    if (op->hasOption("precon_type"))
        _option.precon_type = _option.getPreconType(op->getOption("precon_type"));
    if (op->hasOption("error_tolerance"))
        _option.error_tolerance = op->getOptionAsNum<double>("error_tolerance");
    if (op->hasOption("max_iteration_step"))
        _option.max_iteration_step = op->getOption<int>("max_iteration_step");
}

void SparseLinearEquation::solveEqs(CRSMatrix<double, unsigned> *A, double *rhs, double *x)
{
    double eps = _option.error_tolerance;
    std::size_t steps =  _option.max_iteration_step;
    switch (_option.solver_type)
    {
    case SolverCG:
        CG(A, rhs, x, eps, steps);
        std::cout << "MathLib::CG converged within " << steps << ", residuum is " << eps << std::endl;
        break;
    case SolverBiCGStab:
        BiCGStab(*A, rhs, x, eps, steps);
        break;
    default:
        break;
    }
}

} // end namespace

