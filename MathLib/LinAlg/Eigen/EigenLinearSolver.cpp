/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLinearSolver.h"

#include <boost/property_tree/ptree.hpp>
#include <logog/include/logog.hpp>
#include <unsupported/Eigen/IterativeSolvers>

#include "EigenVector.h"
#include "EigenMatrix.h"
#include "EigenTools.h"

namespace MathLib
{

using boost::property_tree::ptree;

namespace details
{

template <class T_SOLVER>
class EigenDirectSolver : public EigenLinearSolver::IEigenSolver
{
public:
    EigenDirectSolver(EigenMatrix::RawMatrixType &A) : _A(A)
    {
        // fill A and b;
        // Compute the ordering permutation vector from the structural pattern of A
        _solver.analyzePattern(A);
        // Compute the numerical factorization
        _solver.factorize(A);
        if(_solver.info()!=Eigen::Success) {
            ERR("The numerical factorization failed in Eigen");
            return;
        }
    }

    virtual ~EigenDirectSolver() {}

    void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &) override
    {
        //Use the factors to solve the linear system
        INFO("-> solve");
        x = _solver.solve(b);
    }

private:
    T_SOLVER _solver;
    EigenMatrix::RawMatrixType& _A;
};

template <class T_SOLVER>
class EigenIterativeSolver : public EigenLinearSolver::IEigenSolver
{
public:
    EigenIterativeSolver(EigenMatrix::RawMatrixType &A) : _A(A)
    {
        INFO("-> initialize with the coefficient matrix");
        _solver.compute(A);
        if(_solver.info()!=Eigen::Success) {
            INFO("\t failed");
            return;
        }
    }

    virtual ~EigenIterativeSolver() {}

    void solve(EigenVector::RawVectorType &b, EigenVector::RawVectorType &x, EigenOption &opt) override
    {
        INFO("-> solve");
        _solver.setTolerance(opt.error_tolerance);
        _solver.setMaxIterations(opt.max_iterations);
        x = _solver.solveWithGuess(b, x);
        if(_solver.info()!=Eigen::Success) {
            INFO("\t solving failed");
          return;
        }
        INFO("\t iteration: %d/%ld", _solver.iterations(), opt.max_iterations);
        INFO("\t residual: %e\n", _solver.error());
    }

private:
    T_SOLVER _solver;
    EigenMatrix::RawMatrixType& _A;
};

} // details

EigenLinearSolver::EigenLinearSolver(EigenMatrix &A, ptree const*const option)
{
    if (option)
        setOption(*option);

    A.getRawMatrix().makeCompressed();
    if (_option.solver_type==EigenOption::SolverType::SparseLU) {
        typedef Eigen::SparseLU<EigenMatrix::RawMatrixType, Eigen::COLAMDOrdering<int> > SolverType;
        _solver = new details::EigenDirectSolver<SolverType>(A.getRawMatrix());
    } else if (_option.solver_type==EigenOption::SolverType::BiCGSTAB) {
        typedef Eigen::BiCGSTAB<EigenMatrix::RawMatrixType, Eigen::DiagonalPreconditioner<double>> SolverType;
        _solver = new details::EigenIterativeSolver<SolverType>(A.getRawMatrix());
    } else if (_option.solver_type==EigenOption::SolverType::CG) {
        typedef Eigen::ConjugateGradient<EigenMatrix::RawMatrixType, Eigen::Lower, Eigen::DiagonalPreconditioner<double>> SolverType;
        _solver = new details::EigenIterativeSolver<SolverType>(A.getRawMatrix());
    }
}

void EigenLinearSolver::setOption(const ptree &option)
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

void EigenLinearSolver::solve(EigenVector &b, EigenVector &x)
{
    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver computation");
    _solver->solve(b.getRawVector(), x.getRawVector(), _option);
    INFO("------------------------------------------------------------------");
}

} //MathLib
