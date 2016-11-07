/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLinearSolver.h"

#include <logog/include/logog.hpp>

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#include "BaseLib/ConfigTree.h"
#include "EigenVector.h"
#include "EigenMatrix.h"
#include "EigenTools.h"

#include "MathLib/LinAlg/LinearSolverOptions.h"

namespace MathLib
{

// TODO change to LinearSolver
class EigenLinearSolverBase
{
public:
    using Vector = EigenVector::RawVectorType;
    using Matrix = EigenMatrix::RawMatrixType;

    virtual ~EigenLinearSolverBase() = default;

    //! Solves the linear equation system \f$ A x = b \f$ for \f$ x \f$.
    virtual bool solve(Matrix &A, Vector const& b, Vector &x, EigenOption &opt) = 0;
};

namespace details
{

/// Template class for Eigen direct linear solvers
template <class T_SOLVER>
class EigenDirectLinearSolver final : public EigenLinearSolverBase
{
public:
    bool solve(Matrix& A, Vector const& b, Vector& x, EigenOption& /*opt*/) override
    {
        INFO("-> solve");
        if (!A.isCompressed()) A.makeCompressed();

        _solver.compute(A);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solver initialization");
            return false;
        }

        x = _solver.solve(b);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        return true;
    }

private:
    T_SOLVER _solver;
};

/// Template class for Eigen iterative linear solvers
template <class T_SOLVER>
class EigenIterativeLinearSolver final : public EigenLinearSolverBase
{
public:
    bool solve(Matrix& A, Vector const& b, Vector& x, EigenOption& opt) override
    {
        INFO("-> solve");
        _solver.setTolerance(opt.error_tolerance);
        _solver.setMaxIterations(opt.max_iterations);

        if (!A.isCompressed())
            A.makeCompressed();

        _solver.compute(A);
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solver initialization");
            return false;
        }

        x = _solver.solveWithGuess(b, x);
        INFO("\t iteration: %d/%ld", _solver.iterations(), opt.max_iterations);
        INFO("\t residual: %e\n", _solver.error());

        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        return true;
    }

private:
    T_SOLVER _solver;
};

template <template <typename, typename> class Solver, typename Precon>
std::unique_ptr<EigenLinearSolverBase> createIterativeSolver()
{
    using Slv = EigenIterativeLinearSolver<
        Solver<EigenMatrix::RawMatrixType, Precon>>;
    return std::unique_ptr<Slv>(new Slv);
}

template <template <typename, typename> class Solver>
std::unique_ptr<EigenLinearSolverBase> createIterativeSolver(
    EigenOption::PreconType precon_type)
{
    switch (precon_type) {
        case EigenOption::PreconType::NONE:
            return createIterativeSolver<Solver,
                                         Eigen::IdentityPreconditioner>();
        case EigenOption::PreconType::DIAGONAL:
            return createIterativeSolver<
                Solver, Eigen::DiagonalPreconditioner<double>>();
        case EigenOption::PreconType::ILUT:
            // TODO for this preconditioner further options can be passed.
            // see https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html
            return createIterativeSolver<
                Solver, Eigen::IncompleteLUT<double>>();
        default:
            OGS_FATAL("Invalid Eigen preconditioner type.");
    }
}

template <typename Mat, typename Precon>
using EigenCGSolver = Eigen::ConjugateGradient<Mat, Eigen::Lower, Precon>;

std::unique_ptr<EigenLinearSolverBase> createIterativeSolver(
    EigenOption::SolverType solver_type, EigenOption::PreconType precon_type)
{
    switch (solver_type) {
        case EigenOption::SolverType::BiCGSTAB: {
            return createIterativeSolver<Eigen::BiCGSTAB>(precon_type);
        }
        case EigenOption::SolverType::CG: {
            return createIterativeSolver<EigenCGSolver>(precon_type);
        }
        default:
            OGS_FATAL("Invalid Eigen iterative linear solver type. Aborting.");
    }
}

} // details

EigenLinearSolver::EigenLinearSolver(
                            const std::string& /*solver_name*/,
                            const BaseLib::ConfigTree* const option)
{
    using Matrix = EigenMatrix::RawMatrixType;

    if (option)
        setOption(*option);

    // TODO for my taste it is much too unobvious that the default solver type
    //      currently is SparseLU.
    switch (_option.solver_type) {
        case EigenOption::SolverType::SparseLU: {
            using SolverType =
                Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>>;
            _solver.reset(new details::EigenDirectLinearSolver<SolverType>);
            return;
        }
        case EigenOption::SolverType::BiCGSTAB:
        case EigenOption::SolverType::CG:
            _solver = details::createIterativeSolver(_option.solver_type,
                                                     _option.precon_type);
            return;
        case EigenOption::SolverType::PardisoLU: {
#ifdef USE_MKL
            using SolverType = Eigen::PardisoLU<EigenMatrix::RawMatrixType>;
            _solver.reset(new details::EigenDirectLinearSolver<SolverType>);
            return;
#else
            OGS_FATAL(
                "The code is not compiled with Intel MKL. Linear solver type "
                "PardisoLU is not available.");
#endif
        }
    }

    OGS_FATAL("Invalid Eigen linear solver type. Aborting.");
}

EigenLinearSolver::~EigenLinearSolver() = default;

void EigenLinearSolver::setOption(BaseLib::ConfigTree const& option)
{
    ignoreOtherLinearSolvers(option, "eigen");
    //! \ogs_file_param{linear_solver__eigen}
    auto const ptSolver = option.getConfigSubtreeOptional("eigen");
    if (!ptSolver)
        return;

    if (auto solver_type =
            //! \ogs_file_param{linear_solver__eigen__solver_type}
            ptSolver->getConfigParameterOptional<std::string>("solver_type")) {
        _option.solver_type = _option.getSolverType(*solver_type);
    }
    if (auto precon_type =
            //! \ogs_file_param{linear_solver__eigen__precon_type}
            ptSolver->getConfigParameterOptional<std::string>("precon_type")) {
        _option.precon_type = _option.getPreconType(*precon_type);
    }
    if (auto error_tolerance =
            //! \ogs_file_param{linear_solver__eigen__error_tolerance}
            ptSolver->getConfigParameterOptional<double>("error_tolerance")) {
        _option.error_tolerance = *error_tolerance;
    }
    if (auto max_iteration_step =
            //! \ogs_file_param{linear_solver__eigen__max_iteration_step}
            ptSolver->getConfigParameterOptional<int>("max_iteration_step")) {
        _option.max_iterations = *max_iteration_step;
    }
}

bool EigenLinearSolver::solve(EigenMatrix &A, EigenVector& b, EigenVector &x)
{
    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver computation");

    auto const success = _solver->solve(A.getRawMatrix(), b.getRawVector(),
                                        x.getRawVector(), _option);

    INFO("------------------------------------------------------------------");

    return success;
}

} //MathLib
