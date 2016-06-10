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
        if(_solver.info()!=Eigen::Success) {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        INFO("\t iteration: %d/%ld", _solver.iterations(), opt.max_iterations);
        INFO("\t residual: %e\n", _solver.error());

        return true;
    }

private:
    T_SOLVER _solver;
};

} // details

EigenLinearSolver::EigenLinearSolver(
                            const std::string /*solver_name*/,
                            const BaseLib::ConfigTree* const option)
{
    using Matrix = EigenMatrix::RawMatrixType;

    if (option)
        setOption(*option);

    // TODO for my taste it is much too unobvious that the default solver type
    //      currently is SparseLU.
    switch (_option.solver_type)
    {
    case EigenOption::SolverType::SparseLU: {
        using SolverType = Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>>;
        _solver.reset(new details::EigenDirectLinearSolver<SolverType>);
        break;
    }
    case EigenOption::SolverType::BiCGSTAB: {
        using SolverType = Eigen::BiCGSTAB<Matrix, Eigen::DiagonalPreconditioner<double>>;
        _solver.reset(new details::EigenIterativeLinearSolver<SolverType>);
        break;
    }
    case EigenOption::SolverType::CG: {
        using SolverType = Eigen::ConjugateGradient<Matrix, Eigen::Lower, Eigen::DiagonalPreconditioner<double>>;
        _solver.reset(new details::EigenIterativeLinearSolver<SolverType>);
        break;
    }
    case EigenOption::SolverType::INVALID:
        ERR("Invalid Eigen linear solver type. Aborting.");
        std::abort();
    }
}

EigenLinearSolver::~EigenLinearSolver() = default;

void EigenLinearSolver::setOption(BaseLib::ConfigTree const& option)
{
    ignoreOtherLinearSolvers(option, "eigen");
    //! \ogs_file_param{linear_solver__eigen}
    auto const ptSolver = option.getSubtreeOptional("eigen");
    if (!ptSolver)
        return;

    //! \ogs_file_param{linear_solver__eigen__solver_type}
    if (auto solver_type = ptSolver->getParameterOptional<std::string>("solver_type")) {
        _option.solver_type = _option.getSolverType(*solver_type);
    }
    //! \ogs_file_param{linear_solver__eigen__precon_type}
    if (auto precon_type = ptSolver->getParameterOptional<std::string>("precon_type")) {
        _option.precon_type = _option.getPreconType(*precon_type);
    }
    //! \ogs_file_param{linear_solver__eigen__error_tolerance}
    if (auto error_tolerance = ptSolver->getParameterOptional<double>("error_tolerance")) {
        _option.error_tolerance = *error_tolerance;
    }
    //! \ogs_file_param{linear_solver__eigen__max_iteration_step}
    if (auto max_iteration_step = ptSolver->getParameterOptional<int>("max_iteration_step")) {
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
