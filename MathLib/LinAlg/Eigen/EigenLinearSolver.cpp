/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLinearSolver.h"

#include <Eigen/Sparse>

#include "BaseLib/Logging.h"

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#ifdef USE_EIGEN_UNSUPPORTED
#include <unsupported/Eigen/src/IterativeSolvers/GMRES.h>
#include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>
#endif

#include "BaseLib/ConfigTree.h"
#include "EigenMatrix.h"
#include "EigenTools.h"
#include "EigenVector.h"
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
    virtual bool solve(Matrix& A, Vector const& b, Vector& x,
                       EigenOption& opt) = 0;
};

namespace details
{
/// Template class for Eigen direct linear solvers
template <class T_SOLVER>
class EigenDirectLinearSolver final : public EigenLinearSolverBase
{
public:
    bool solve(Matrix& A, Vector const& b, Vector& x, EigenOption& opt) override
    {
        INFO("-> solve with {:s}", EigenOption::getSolverName(opt.solver_type));
        if (!A.isCompressed())
        {
            A.makeCompressed();
        }

        solver_.compute(A);
        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solver initialization");
            return false;
        }

        x = solver_.solve(b);
        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        return true;
    }

private:
    T_SOLVER solver_;
};

/// Template class for Eigen iterative linear solvers
template <class T_SOLVER>
class EigenIterativeLinearSolver final : public EigenLinearSolverBase
{
public:
    bool solve(Matrix& A, Vector const& b, Vector& x, EigenOption& opt) override
    {
        INFO("-> solve with {:s} (precon {:s})",
             EigenOption::getSolverName(opt.solver_type),
             EigenOption::getPreconName(opt.precon_type));
        solver_.setTolerance(opt.error_tolerance);
        solver_.setMaxIterations(opt.max_iterations);
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setRestart(
            opt.restart);

        if (!A.isCompressed())
        {
            A.makeCompressed();
        }

        solver_.compute(A);
        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solver initialization");
            return false;
        }

        x = solver_.solveWithGuess(b, x);
        INFO("\t iteration: {:d}/{:d}", solver_.iterations(),
             opt.max_iterations);
        INFO("\t residual: {:e}\n", solver_.error());

        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        return true;
    }

private:
    T_SOLVER solver_;
    void setRestart(int const /*restart*/) {}
};

/// Specialization for (all) three preconditioners separately
template <>
void EigenIterativeLinearSolver<
    Eigen::GMRES<EigenMatrix::RawMatrixType,
                 Eigen::IdentityPreconditioner>>::setRestart(int const restart)
{
    solver_.set_restart(restart);
    INFO("-> set restart value: {:d}", solver_.get_restart());
}

template <>
void EigenIterativeLinearSolver<Eigen::GMRES<
    EigenMatrix::RawMatrixType,
    Eigen::DiagonalPreconditioner<double>>>::setRestart(int const restart)
{
    solver_.set_restart(restart);
    INFO("-> set restart value: {:d}", solver_.get_restart());
}

template <>
void EigenIterativeLinearSolver<
    Eigen::GMRES<EigenMatrix::RawMatrixType,
                 Eigen::IncompleteLUT<double>>>::setRestart(int const restart)
{
    solver_.set_restart(restart);
    INFO("-> set restart value: {:d}", solver_.get_restart());
}

template <template <typename, typename> class Solver, typename Precon>
std::unique_ptr<EigenLinearSolverBase> createIterativeSolver()
{
    using Slv =
        EigenIterativeLinearSolver<Solver<EigenMatrix::RawMatrixType, Precon>>;
    return std::make_unique<Slv>();
}

template <template <typename, typename> class Solver>
std::unique_ptr<EigenLinearSolverBase> createIterativeSolver(
    EigenOption::PreconType precon_type)
{
    switch (precon_type)
    {
        case EigenOption::PreconType::NONE:
            return createIterativeSolver<Solver,
                                         Eigen::IdentityPreconditioner>();
        case EigenOption::PreconType::DIAGONAL:
            return createIterativeSolver<
                Solver, Eigen::DiagonalPreconditioner<double>>();
        case EigenOption::PreconType::ILUT:
            // TODO for this preconditioner further options can be passed.
            // see
            // https://eigen.tuxfamily.org/dox/classEigen_1_1IncompleteLUT.html
            return createIterativeSolver<Solver,
                                         Eigen::IncompleteLUT<double>>();
        default:
            OGS_FATAL("Invalid Eigen preconditioner type.");
    }
}

template <typename Mat, typename Precon>
using EigenCGSolver = Eigen::ConjugateGradient<Mat, Eigen::Lower, Precon>;

std::unique_ptr<EigenLinearSolverBase> createIterativeSolver(
    EigenOption::SolverType solver_type, EigenOption::PreconType precon_type)
{
    switch (solver_type)
    {
        case EigenOption::SolverType::BiCGSTAB:
        {
            return createIterativeSolver<Eigen::BiCGSTAB>(precon_type);
        }
        case EigenOption::SolverType::CG:
        {
            return createIterativeSolver<EigenCGSolver>(precon_type);
        }
        case EigenOption::SolverType::GMRES:
        {
#ifdef USE_EIGEN_UNSUPPORTED
            return createIterativeSolver<Eigen::GMRES>(precon_type);
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "Linear solver type GMRES is not available.");
#endif
        }
        default:
            OGS_FATAL("Invalid Eigen iterative linear solver type. Aborting.");
    }
}

}  // namespace details

EigenLinearSolver::EigenLinearSolver(const std::string& /*solver_name*/,
                                     const BaseLib::ConfigTree* const option)
{
    using Matrix = EigenMatrix::RawMatrixType;

    if (option)
    {
        setOption(*option);
    }

    // TODO for my taste it is much too unobvious that the default solver type
    //      currently is SparseLU.
    switch (option_.solver_type)
    {
        case EigenOption::SolverType::SparseLU:
        {
            using SolverType =
                Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>>;
            solver_ = std::make_unique<
                details::EigenDirectLinearSolver<SolverType>>();
            return;
        }
        case EigenOption::SolverType::BiCGSTAB:
        case EigenOption::SolverType::CG:
        case EigenOption::SolverType::GMRES:
            solver_ = details::createIterativeSolver(option_.solver_type,
                                                     option_.precon_type);
            return;
        case EigenOption::SolverType::PardisoLU:
        {
#ifdef USE_MKL
            using SolverType = Eigen::PardisoLU<EigenMatrix::RawMatrixType>;
            solver_.reset(new details::EigenDirectLinearSolver<SolverType>);
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

EigenLinearSolver::EigenLinearSolver(std::string const& /*solver_name*/,
                                     EigenOption const& option)
    : option_(option)
{
    using Matrix = EigenMatrix::RawMatrixType;

    switch (option_.solver_type)
    {
        case EigenOption::SolverType::SparseLU:
        {
            using SolverType =
                Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int>>;
            solver_ = std::make_unique<
                details::EigenDirectLinearSolver<SolverType>>();
            return;
        }
        case EigenOption::SolverType::BiCGSTAB:
        case EigenOption::SolverType::CG:
        case EigenOption::SolverType::GMRES:
            solver_ = details::createIterativeSolver(option_.solver_type,
                                                     option_.precon_type);
            return;
        case EigenOption::SolverType::PardisoLU:
        {
#ifdef USE_MKL
            using SolverType = Eigen::PardisoLU<EigenMatrix::RawMatrixType>;
            solver_.reset(new details::EigenDirectLinearSolver<SolverType>);
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
    //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen}
    auto const ptSolver = option.getConfigSubtreeOptional("eigen");
    if (!ptSolver)
    {
        return;
    }

    if (auto solver_type =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__solver_type}
        ptSolver->getConfigParameterOptional<std::string>("solver_type"))
    {
        option_.solver_type = MathLib::EigenOption::getSolverType(*solver_type);
    }
    if (auto precon_type =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__precon_type}
        ptSolver->getConfigParameterOptional<std::string>("precon_type"))
    {
        option_.precon_type = MathLib::EigenOption::getPreconType(*precon_type);
    }
    if (auto error_tolerance =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__error_tolerance}
        ptSolver->getConfigParameterOptional<double>("error_tolerance"))
    {
        option_.error_tolerance = *error_tolerance;
    }
    if (auto max_iteration_step =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__max_iteration_step}
        ptSolver->getConfigParameterOptional<int>("max_iteration_step"))
    {
        option_.max_iterations = *max_iteration_step;
    }
    if (auto scaling =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__scaling}
        ptSolver->getConfigParameterOptional<bool>("scaling"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        option_.scaling = *scaling;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules. "
            "scaling is not available.");
#endif
    }
    if (auto restart =
            //! \ogs_file_param{prj__linear_solvers__linear_solver__eigen__restart}
        ptSolver->getConfigParameterOptional<int>("restart"))
    {
#ifdef USE_EIGEN_UNSUPPORTED
        option_.restart = *restart;
#else
        OGS_FATAL(
            "The code is not compiled with the Eigen unsupported modules. "
            "GMRES/GMRES option restart is not available.");
#endif
    }
}

bool EigenLinearSolver::solve(EigenMatrix& A, EigenVector& b, EigenVector& x)
{
    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver computation");

#ifdef USE_EIGEN_UNSUPPORTED
    std::unique_ptr<Eigen::IterScaling<EigenMatrix::RawMatrixType>> scal;
    if (option_.scaling)
    {
        INFO("-> scale");
        scal =
            std::make_unique<Eigen::IterScaling<EigenMatrix::RawMatrixType>>();
        scal->computeRef(A.getRawMatrix());
        b.getRawVector() = scal->LeftScaling().cwiseProduct(b.getRawVector());
    }
#endif
    auto const success = solver_->solve(A.getRawMatrix(), b.getRawVector(),
                                        x.getRawVector(), option_);
#ifdef USE_EIGEN_UNSUPPORTED
    if (scal)
    {
        x.getRawVector() = scal->RightScaling().cwiseProduct(x.getRawVector());
    }
#endif

    INFO("------------------------------------------------------------------");

    return success;
}

}  // namespace MathLib
