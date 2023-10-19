/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "EigenLinearSolver.h"

#include <Eigen/Sparse>

#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

#ifdef USE_EIGEN_UNSUPPORTED
// TODO(naumov); Don't include header files directly from Eigen's source dir
// clang-format off
#include <unsupported/Eigen/IterativeSolvers>
#include <unsupported/Eigen/src/IterativeSolvers/Scaling.h>
// clang-format on
#endif

#include "EigenMatrix.h"
#include "EigenVector.h"

namespace MathLib
{
// TODO change to LinearSolver
class EigenLinearSolverBase
{
public:
    using Vector = EigenVector::RawVectorType;
    using Matrix = EigenMatrix::RawMatrixType;

    virtual ~EigenLinearSolverBase() = default;

    bool solve(Vector& b, Vector& x, EigenOption& opt)
    {
#ifdef USE_EIGEN_UNSUPPORTED
        if (scaling_)
        {
            b = scaling_->LeftScaling().cwiseProduct(b);
        }
#endif

        auto const success = solveImpl(b, x, opt);

#ifdef USE_EIGEN_UNSUPPORTED
        if (scaling_)
        {
            x = scaling_->RightScaling().cwiseProduct(x);
        }
#endif

        return success;
    }

    bool compute(Matrix& A, EigenOption& opt)
    {
#ifdef USE_EIGEN_UNSUPPORTED
        if (opt.scaling)
        {
            INFO("-> scale");
            scaling_ = std::make_unique<
                Eigen::IterScaling<EigenMatrix::RawMatrixType>>();
            scaling_->computeRef(A);
        }
#endif

        return computeImpl(A, opt);
    }

protected:
    virtual bool solveImpl(Vector const& b, Vector& x, EigenOption& opt) = 0;

    virtual bool computeImpl(Matrix& A, EigenOption& opt) = 0;

private:
#ifdef USE_EIGEN_UNSUPPORTED
    std::unique_ptr<Eigen::IterScaling<EigenMatrix::RawMatrixType>> scaling_;
#endif
};

namespace details
{
/// Template class for Eigen direct linear solvers
template <class T_SOLVER>
class EigenDirectLinearSolver final : public EigenLinearSolverBase
{
public:
    bool solveImpl(Vector const& b, Vector& x, EigenOption& opt) override
    {
        INFO("-> solve with Eigen direct linear solver {:s}",
             EigenOption::getSolverName(opt.solver_type));

        x = solver_.solve(b);
        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solve");
            return false;
        }

        return true;
    }

    bool computeImpl(Matrix& A, EigenOption& opt) override
    {
        INFO("-> compute with Eigen direct linear solver {:s}",
             EigenOption::getSolverName(opt.solver_type));
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

        return true;
    }

private:
    T_SOLVER solver_;
};

#ifdef USE_EIGEN_UNSUPPORTED
// implementations for some iterative linear solver methods --------------------

// restart
template <typename Solver>
void setRestartImpl(Solver&, int const)
{
    DBUG("-> restart is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setRestartImpl(Eigen::GMRES<Matrix, Precon>& solver, int const restart)
{
    solver.set_restart(restart);
    INFO("-> set restart value: {:d}", solver.get_restart());
}

// L
template <typename Solver>
void setLImpl(Solver&, int const)
{
    DBUG("-> setL() is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setLImpl(Eigen::BiCGSTABL<Matrix, Precon>& solver, int const l)
{
    solver.setL(l);
}

template <typename Matrix, typename Precon>
void setLImpl(Eigen::IDRSTABL<Matrix, Precon>& solver, int const l)
{
    solver.setL(l);
}

// S
template <typename Solver>
void setSImpl(Solver&, int const)
{
    DBUG("-> setS() is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setSImpl(Eigen::IDRS<Matrix, Precon>& solver, int const s)
{
    solver.setS(s);
}

template <typename Matrix, typename Precon>
void setSImpl(Eigen::IDRSTABL<Matrix, Precon>& solver, int const s)
{
    solver.setS(s);
}

// angle
template <typename Solver>
void setAngleImpl(Solver&, double const)
{
    DBUG("-> setAngle() is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setAngleImpl(Eigen::IDRS<Matrix, Precon>& solver, double const angle)
{
    solver.setAngle(angle);
}

// smoothing
template <typename Solver>
void setSmoothingImpl(Solver&, bool const)
{
    DBUG("-> setSmoothing() is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setSmoothingImpl(Eigen::IDRS<Matrix, Precon>& solver, bool const smoothing)
{
    solver.setSmoothing(smoothing);
}

// residual update
template <typename Solver>
void setResidualUpdateImpl(Solver&, bool const)
{
    DBUG("-> setResidualUpdate() is not implemented for this linear solver.");
}

template <typename Matrix, typename Precon>
void setResidualUpdateImpl(Eigen::IDRS<Matrix, Precon>& solver,
                           bool const residual_update)
{
    solver.setResidualUpdate(residual_update);
}

// -----------------------------------------------------------------------------
#endif

/// Template class for Eigen iterative linear solvers
template <class T_SOLVER>
class EigenIterativeLinearSolver final : public EigenLinearSolverBase
{
public:
    bool computeImpl(Matrix& A, EigenOption& opt) override
    {
        INFO("-> compute with Eigen iterative linear solver {:s} (precon {:s})",
             EigenOption::getSolverName(opt.solver_type),
             EigenOption::getPreconName(opt.precon_type));
        solver_.setTolerance(opt.error_tolerance);
        solver_.setMaxIterations(opt.max_iterations);

#ifdef USE_EIGEN_UNSUPPORTED
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setRestart(
            opt.restart);
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setL(opt.l);
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setS(opt.s);
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setSmoothing(
            opt.smoothing);
        MathLib::details::EigenIterativeLinearSolver<T_SOLVER>::setAngle(
            opt.angle);
        MathLib::details::EigenIterativeLinearSolver<
            T_SOLVER>::setResidualUpdate(opt.residualupdate);
#endif

        // matrix must be copied, because Eigen's linear solver stores a
        // reference to it cf.
        // https://libeigen.gitlab.io/docs/classEigen_1_1IterativeSolverBase.html#a7dfa55c55e82d697bde227696a630914
        A_ = A;

        if (!A_.isCompressed())
        {
            A_.makeCompressed();
        }

        solver_.compute(A_);
        if (solver_.info() != Eigen::Success)
        {
            ERR("Failed during Eigen linear solver initialization");
            return false;
        }

        return true;
    }

    bool solveImpl(Vector const& b, Vector& x, EigenOption& opt) override
    {
        INFO("-> solve with Eigen iterative linear solver {:s} (precon {:s})",
             EigenOption::getSolverName(opt.solver_type),
             EigenOption::getPreconName(opt.precon_type));

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
    Matrix A_;
#ifdef USE_EIGEN_UNSUPPORTED
    void setRestart(int const restart) { setRestartImpl(solver_, restart); }
    void setL(int const l) { setLImpl(solver_, l); }
    void setS(int const s) { setSImpl(solver_, s); }
    void setAngle(double const angle) { setAngleImpl(solver_, angle); }
    void setSmoothing(bool const smoothing)
    {
        setSmoothingImpl(solver_, smoothing);
    }
    void setResidualUpdate(bool const residual_update)
    {
        setResidualUpdateImpl(solver_, residual_update);
    }
#endif
};

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
            // https://libeigen.gitlab.io/docs/classEigen_1_1IncompleteLUT.html
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
        case EigenOption::SolverType::BiCGSTABL:
        {
#ifdef USE_EIGEN_UNSUPPORTED
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
            return createIterativeSolver<Eigen::BiCGSTABL>(precon_type);
#else
            OGS_FATAL("BiCGSTABL requires at least Eigen version 3.4.90");
#endif
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "Linear solver type BiCGSTABL is not available.");
#endif
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
        case EigenOption::SolverType::IDRS:
        {
#ifdef USE_EIGEN_UNSUPPORTED
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
            return createIterativeSolver<Eigen::IDRS>(precon_type);
#else
            OGS_FATAL("IDRS requires at least Eigen version 3.4.90");
#endif
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "Linear solver type IDRS is not available.");
#endif
        }
        case EigenOption::SolverType::IDRSTABL:
        {
#ifdef USE_EIGEN_UNSUPPORTED
#if EIGEN_VERSION_AT_LEAST(3, 4, 90)
            return createIterativeSolver<Eigen::IDRSTABL>(precon_type);
#else
            OGS_FATAL("IDRSTABL requires at least Eigen version 3.4.90");
#endif
#else
            OGS_FATAL(
                "The code is not compiled with the Eigen unsupported modules. "
                "Linear solver type IDRSTABL is not available.");
#endif
        }
        default:
            OGS_FATAL("Invalid Eigen iterative linear solver type. Aborting.");
    }
}

}  // namespace details

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
        case EigenOption::SolverType::BiCGSTABL:
        case EigenOption::SolverType::CG:
        case EigenOption::SolverType::GMRES:
        case EigenOption::SolverType::IDRS:
        case EigenOption::SolverType::IDRSTABL:
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

bool EigenLinearSolver::compute(EigenMatrix& A)
{
    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver compute()");

    return solver_->compute(A.getRawMatrix(), option_);
}

bool EigenLinearSolver::solve(EigenVector& b, EigenVector& x)
{
    INFO("------------------------------------------------------------------");
    INFO("*** Eigen solver solve()");

    return solver_->solve(b.getRawVector(), x.getRawVector(), option_);
}

bool EigenLinearSolver::solve(EigenMatrix& A, EigenVector& b, EigenVector& x)
{
    return solver_->compute(A.getRawMatrix(), option_) &&
           solver_->solve(b.getRawVector(), x.getRawVector(), option_);
}

}  // namespace MathLib
