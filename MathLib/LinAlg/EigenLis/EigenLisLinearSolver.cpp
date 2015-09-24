/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLisLinearSolver.h"

#include <boost/property_tree/ptree.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <logog/include/logog.hpp>

#include "MathLib/LinAlg/Eigen/EigenMatrix.h"
#include "MathLib/LinAlg/Eigen/EigenVector.h"
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#include "MathLib/LinAlg/Lis/LisVector.h"
#include "MathLib/LinAlg/Lis/LisLinearSolver.h"

namespace MathLib
{

using boost::property_tree::ptree;

EigenLisLinearSolver::EigenLisLinearSolver(EigenMatrix &A,
                      const std::string /*solver_name*/,
                      ptree const*const option)
: _A(A)
{
    if (option)
        setOption(*option);
}

void EigenLisLinearSolver::setOption(const ptree &option)
{
    boost::optional<ptree> ptSolver = option.get_child("LinearSolver");
    if (!ptSolver)
        return;

    boost::optional<std::string> solver_type = ptSolver->get_optional<std::string>("solver_type");
    if (solver_type) {
        _option.solver_type = _option.getSolverType(*solver_type);
    }
    boost::optional<std::string> precon_type = ptSolver->get_optional<std::string>("precon_type");
    if (precon_type) {
        _option.precon_type = _option.getPreconType(*precon_type);
    }
    boost::optional<double> error_tolerance = ptSolver->get_optional<double>("error_tolerance");
    if (error_tolerance) {
        _option.error_tolerance = *error_tolerance;
    }
    boost::optional<int> max_iteration_step = ptSolver->get_optional<int>("max_iteration_step");
    if (max_iteration_step) {
        _option.max_iterations = *max_iteration_step;
    }
}

void EigenLisLinearSolver::solve(EigenVector &b_, EigenVector &x_)
{
    auto &A = _A.getRawMatrix();
    auto &b = b_.getRawVector();
    auto &x = x_.getRawVector();

    if (!A.isCompressed())
        A.makeCompressed();
    int nnz = A.nonZeros();
    int* ptr = A.outerIndexPtr();
    int* col = A.innerIndexPtr();
    double* data = A.valuePtr();
    LisMatrix lisA(_A.getNRows(), nnz, ptr, col, data);
    LisVector lisb(b.rows(), b.data());
    LisVector lisx(x.rows(), x.data());

    LisLinearSolver lissol(lisA);
    lissol.setOption(_option);
    lissol.solve(lisb, lisx);

    for (std::size_t i=0; i<lisx.size(); i++)
        x[i] = lisx[i];
}

} //MathLib
