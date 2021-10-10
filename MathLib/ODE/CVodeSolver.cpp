/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "CVodeSolver.h"

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include <cassert>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

//! \addtogroup ExternalODESolverInterface
//! @{

/*! Checks Sundials error flags and aborts the program in case of error.
 *
 * \param f_name name of the function that returned the \c error_flag
 * \param error_flag the error flag to be checked
 */
void check_error(std::string const& f_name, int const error_flag)
{
    if (error_flag != CV_SUCCESS)
    {
        OGS_FATAL("CVodeSolver: {:s} failed with error flag {:d}.", f_name,
                  error_flag);
    }
}

//! Prints some statistics about an ODE solver run.
void printStats(void* cvode_mem)
{
    long int nst = 0, nfe = 0, nsetups = 0, nje = 0, nfeLS = 0, nni = 0,
             ncfn = 0, netf = 0, nge = 0;

    check_error("CVodeGetNumSteps", CVodeGetNumSteps(cvode_mem, &nst));
    check_error("CVodeGetNumRhsEvals", CVodeGetNumRhsEvals(cvode_mem, &nfe));
    check_error("CVodeGetNumLinSolvSetups",
                CVodeGetNumLinSolvSetups(cvode_mem, &nsetups));
    check_error("CVodeGetNumErrTestFails",
                CVodeGetNumErrTestFails(cvode_mem, &netf));
    check_error("CVodeGetNumNonlinSolvIters",
                CVodeGetNumNonlinSolvIters(cvode_mem, &nni));
    check_error("CVodeGetNumNonlinSolvConvFails",
                CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn));
    check_error("CVDlsGetNumJacEvals", CVDlsGetNumJacEvals(cvode_mem, &nje));
    check_error("CVDlsGetNumRhsEvals", CVDlsGetNumRhsEvals(cvode_mem, &nfeLS));
    check_error("CVodeGetNumGEvals", CVodeGetNumGEvals(cvode_mem, &nge));

    DBUG("Sundials CVode solver. Statistics:");
    DBUG(
        "nst = {:<0d}  nfe = {:<0d} nsetups = {:<0d} nfeLS = {:<0d} nje = {:d}",
        nst, nfe, nsetups, nfeLS, nje);
    DBUG("nni = {:<0d} ncfn = {:<0d}    netf = {:<0d}   nge = {:d}\n", nni,
         ncfn, netf, nge);
}

//! @}

namespace MathLib
{
namespace ODE
{
//! \addtogroup ExternalODESolverInterface
//! @{

/*! This class provides concrete access to Sundials' CVode solver.
 *
 * This class is the implementation part in the pimpl idiom used by the
 * CVodeSolver class. Therefore all if this classes methods are only forwarded
 * from CVodeSolver.
 */
class CVodeSolverImpl final
{
    static_assert(std::is_same_v<realtype, double>,
                  "CVode's realtype is not the same as double");

public:
    CVodeSolverImpl(BaseLib::ConfigTree const& config,
                    unsigned const num_equations);

    void setFunction(std::unique_ptr<detail::FunctionHandles>&& f);

    void preSolve();
    bool solve(const double t_end);

    double const* getSolution() const { return NV_DATA_S(y_); }
    double getTime() const { return t_; }
    void getYDot(const double t, double const* const y, double* const y_dot);
    void setTolerance(const double* abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);
    void setIC(const double t0, double const* const y0);

    ~CVodeSolverImpl();

private:
    N_Vector y_ = nullptr;  //!< The solution vector.

    realtype t_;  //! current time

    N_Vector abstol_ = nullptr;  //!< Array of absolute tolerances.
    realtype reltol_;            //!< Relative tolerance

    unsigned num_equations_;  //!< Number of equations in the ODE system.
    void* cvode_mem_;         //!< CVode's internal memory

    //! Function handles that compute \f$\partial \dot y/\partial y\f$
    //! and \f$\dot y\f$.
    std::unique_ptr<detail::FunctionHandles> f_;

    //! The multistep method used for solving the ODE.
    int linear_multistep_method_ = CV_ADAMS;

    //! Either solve via fixed-point iteration or via Newton-Raphson method.
    int nonlinear_solver_iteration_ = CV_FUNCTIONAL;
};

//! @}

CVodeSolverImpl::CVodeSolverImpl(const BaseLib::ConfigTree& config,
                                 const unsigned num_equations)
{
    if (auto const param =
            //! \ogs_file_param{ode_solver__CVODE__linear_multistep_method}
        config.getConfigParameterOptional<std::string>(
            "linear_multistep_method"))
    {
        DBUG("setting linear multistep method (config: {:s})", param->c_str());

        if (*param == "Adams")
        {
            linear_multistep_method_ = CV_ADAMS;
        }
        else if (*param == "BDF")
        {
            linear_multistep_method_ = CV_BDF;
        }
        else
        {
            OGS_FATAL("unknown linear multistep method: {:s}", param->c_str());
        }
    }

    if (auto const param =
            //! \ogs_file_param{ode_solver__CVODE__nonlinear_solver_iteration}
        config.getConfigParameterOptional<std::string>(
            "nonlinear_solver_iteration"))
    {
        DBUG("setting nonlinear solver iteration (config: {:s})",
             param->c_str());

        if (*param == "Functional")
        {
            nonlinear_solver_iteration_ = CV_FUNCTIONAL;
        }
        else if (*param == "Newton")
        {
            nonlinear_solver_iteration_ = CV_NEWTON;
        }
        else
        {
            OGS_FATAL("unknown nonlinear solver iteration: {:s}",
                      param->c_str());
        }
    }

    y_ = N_VNew_Serial(num_equations);
    abstol_ = N_VNew_Serial(num_equations);
    num_equations_ = num_equations;

    cvode_mem_ =
        CVodeCreate(linear_multistep_method_, nonlinear_solver_iteration_);

    if (cvode_mem_ == nullptr || y_ == nullptr || abstol_ == nullptr)
    {
        OGS_FATAL("couldn't allocate storage for CVode solver.");
    }

    auto f_wrapped = [](const realtype t, const N_Vector y, N_Vector ydot,
                        void* function_handles) -> int
    {
        bool successful =
            static_cast<detail::FunctionHandles*>(function_handles)
                ->call(t, NV_DATA_S(y), NV_DATA_S(ydot));
        return successful ? 0 : 1;
    };

    check_error("CVodeInit", CVodeInit(cvode_mem_, f_wrapped, 0.0, y_));
}

void CVodeSolverImpl::setTolerance(const double* abstol, const double reltol)
{
    for (unsigned i = 0; i < num_equations_; ++i)
    {
        NV_Ith_S(abstol_, i) = abstol[i];
    }

    reltol_ = reltol;
}

void CVodeSolverImpl::setTolerance(const double abstol, const double reltol)
{
    for (unsigned i = 0; i < num_equations_; ++i)
    {
        NV_Ith_S(abstol_, i) = abstol;
    }

    reltol_ = reltol;
}

void CVodeSolverImpl::setFunction(std::unique_ptr<detail::FunctionHandles>&& f)
{
    f_ = std::move(f);
    assert(num_equations_ == f_->getNumberOfEquations());
}

void CVodeSolverImpl::setIC(const double t0, double const* const y0)
{
    for (unsigned i = 0; i < num_equations_; ++i)
    {
        NV_Ith_S(y_, i) = y0[i];
    }

    t_ = t0;
}

void CVodeSolverImpl::preSolve()
{
    assert(f_ != nullptr && "ode function handle was not provided");

    // sets initial conditions
    check_error("CVodeReInit", CVodeReInit(cvode_mem_, t_, y_));

    check_error("CVodeSetUserData",
                CVodeSetUserData(cvode_mem_, static_cast<void*>(f_.get())));

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    check_error("CVodeSVtolerances",
                CVodeSVtolerances(cvode_mem_, reltol_, abstol_));

    /* Call CVDense to specify the CVDENSE dense linear solver */
    check_error("CVDense", CVDense(cvode_mem_, num_equations_));

    if (f_->hasJacobian())
    {
        auto df_wrapped = [](const long N, const realtype t, const N_Vector y,
                             const N_Vector ydot, const DlsMat jac,
                             void* function_handles, N_Vector /*tmp1*/,
                             N_Vector /*tmp2*/, N_Vector /*tmp3*/
                             ) -> int
        {
            (void)N;  // prevent warnings during non-debug build
            auto* fh = static_cast<detail::FunctionHandles*>(function_handles);
            assert(N == fh->getNumberOfEquations());

            // Caution: by calling the DENSE_COL() macro we assume that matrices
            // are stored contiguously in memory!
            // See also the header files sundials_direct.h and cvode_direct.h in
            // the Sundials source code. The comments about the macro DENSE_COL
            // in those files indicate that matrices are stored column-wise.
            bool successful = fh->callJacobian(t, NV_DATA_S(y), NV_DATA_S(ydot),
                                               DENSE_COL(jac, 0));
            return successful ? 0 : 1;
        };

        check_error("CVDlsSetDenseJacFn",
                    CVDlsSetDenseJacFn(cvode_mem_, df_wrapped));
    }
}

bool CVodeSolverImpl::solve(const double t_end)
{
    realtype t_reached;
    check_error("CVode solve",
                CVode(cvode_mem_, t_end, y_, &t_reached, CV_NORMAL));
    t_ = t_reached;

    // check_error asserts that t_end == t_reached and that solving the ODE
    // went fine. Otherwise the program will be aborted. Therefore, we don't
    // have to check manually for errors here and can always safely return true.
    return true;
}

void CVodeSolverImpl::getYDot(const double t, double const* const y,
                              double* const y_dot)
{
    assert(f_ != nullptr);
    f_->call(t, y, y_dot);
}

CVodeSolverImpl::~CVodeSolverImpl()
{
    printStats(cvode_mem_);

    N_VDestroy_Serial(y_);
    N_VDestroy_Serial(abstol_);
    CVodeFree(&cvode_mem_);
}

CVodeSolver::CVodeSolver(BaseLib::ConfigTree const& config,
                         unsigned const num_equations)
    : impl_{new CVodeSolverImpl{config, num_equations}}
{
}

void CVodeSolver::setTolerance(const double* abstol, const double reltol)
{
    impl_->setTolerance(abstol, reltol);
}

void CVodeSolver::setTolerance(const double abstol, const double reltol)
{
    impl_->setTolerance(abstol, reltol);
}

void CVodeSolver::setFunction(std::unique_ptr<detail::FunctionHandles>&& f)
{
    impl_->setFunction(std::move(f));
}

void CVodeSolver::setIC(const double t0, double const* const y0)
{
    impl_->setIC(t0, y0);
}

void CVodeSolver::preSolve()
{
    impl_->preSolve();
}

bool CVodeSolver::solve(const double t_end)
{
    return impl_->solve(t_end);
}

double const* CVodeSolver::getSolution() const
{
    return impl_->getSolution();
}

void CVodeSolver::getYDot(const double t, double const* const y,
                          double* const y_dot) const
{
    impl_->getYDot(t, y, y_dot);
}

double CVodeSolver::getTime() const
{
    return impl_->getTime();
}

CVodeSolver::~CVodeSolver() = default;

}  // namespace ODE
}  // namespace MathLib
