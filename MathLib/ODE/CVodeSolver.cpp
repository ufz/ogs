/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifdef CVODE_FOUND

#include "CVodeSolver.h"

#include <cassert>
#include <logog/include/logog.hpp>

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

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
        OGS_FATAL("CVodeSolver: %s failed with error flag %d.", f_name.c_str(),
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
    DBUG("nst = %-6ld  nfe = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld",
         nst, nfe, nsetups, nfeLS, nje);
    DBUG("nni = %-6ld ncfn = %-6ld    netf = %-6ld   nge = %ld\n", nni, ncfn,
         netf, nge);
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
    static_assert(std::is_same<realtype, double>::value,
                  "CVode's realtype is not the same as double");

public:
    CVodeSolverImpl(BaseLib::ConfigTree const& config,
                    unsigned const num_equations);

    void setFunction(std::unique_ptr<detail::FunctionHandles>&& f);

    void preSolve();
    bool solve(const double t_end);

    double const* getSolution() const { return NV_DATA_S(_y); }
    double getTime() const { return _t; }
    void getYDot(const double t, double const* const y, double* const y_dot);
    void setTolerance(const double* abstol, const double reltol);
    void setTolerance(const double abstol, const double reltol);
    void setIC(const double t0, double const* const y0);

    ~CVodeSolverImpl();

private:
    N_Vector _y = nullptr;  //!< The solution vector.

    realtype _t;  //! current time

    N_Vector _abstol = nullptr;  //!< Array of absolute tolerances.
    realtype _reltol;            //!< Relative tolerance

    unsigned _num_equations;  //!< Number of equations in the ODE system.
    void* _cvode_mem;         //!< CVode's internal memory

    //! Function handles that compute \f$\partial \dot y/\partial y\f$
    //! and \f$\dot y\f$.
    std::unique_ptr<detail::FunctionHandles> _f;

    //! The multistep method used for solving the ODE.
    int _linear_multistep_method = CV_ADAMS;

    //! Either solve via fixed-point iteration or via Newton-Raphson method.
    int _nonlinear_solver_iteration = CV_FUNCTIONAL;
};

//! @}

CVodeSolverImpl::CVodeSolverImpl(const BaseLib::ConfigTree& config,
                                 const unsigned num_equations)
{
    if (auto const param =
        //! \ogs_file_param{ode_solver__CVODE__linear_multistep_method}
        config.getConfigParameterOptional<std::string>("linear_multistep_method"))
    {
        DBUG("setting linear multistep method (config: %s)", param->c_str());

        if (*param == "Adams")
        {
            _linear_multistep_method = CV_ADAMS;
        }
        else if (*param == "BDF")
        {
            _linear_multistep_method = CV_BDF;
        }
        else
        {
            OGS_FATAL("unknown linear multistep method: %s", param->c_str());
        }
    }

    if (auto const param =
        //! \ogs_file_param{ode_solver__CVODE__nonlinear_solver_iteration}
        config.getConfigParameterOptional<std::string>("nonlinear_solver_iteration"))
    {
        DBUG("setting nonlinear solver iteration (config: %s)", param->c_str());

        if (*param == "Functional")
        {
            _nonlinear_solver_iteration = CV_FUNCTIONAL;
        }
        else if (*param == "Newton")
        {
            _nonlinear_solver_iteration = CV_NEWTON;
        }
        else
        {
            OGS_FATAL("unknown nonlinear solver iteration: %s", param->c_str());
        }
    }

    _y = N_VNew_Serial(num_equations);
    _abstol = N_VNew_Serial(num_equations);
    _num_equations = num_equations;

    _cvode_mem =
        CVodeCreate(_linear_multistep_method, _nonlinear_solver_iteration);

    if (_cvode_mem == nullptr || _y == nullptr || _abstol == nullptr)
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

    check_error("CVodeInit", CVodeInit(_cvode_mem, f_wrapped, 0.0, _y));
}

void CVodeSolverImpl::setTolerance(const double* abstol, const double reltol)
{
    for (unsigned i = 0; i < _num_equations; ++i)
    {
        NV_Ith_S(_abstol, i) = abstol[i];
    }

    _reltol = reltol;
}

void CVodeSolverImpl::setTolerance(const double abstol, const double reltol)
{
    for (unsigned i = 0; i < _num_equations; ++i)
    {
        NV_Ith_S(_abstol, i) = abstol;
    }

    _reltol = reltol;
}

void CVodeSolverImpl::setFunction(std::unique_ptr<detail::FunctionHandles>&& f)
{
    _f = std::move(f);
    assert(_num_equations == _f->getNumberOfEquations());
}

void CVodeSolverImpl::setIC(const double t0, double const* const y0)
{
    for (unsigned i = 0; i < _num_equations; ++i)
    {
        NV_Ith_S(_y, i) = y0[i];
    }

    _t = t0;
}

void CVodeSolverImpl::preSolve()
{
    assert(_f != nullptr && "ode function handle was not provided");

    // sets initial conditions
    check_error("CVodeReInit", CVodeReInit(_cvode_mem, _t, _y));

    check_error("CVodeSetUserData",
                CVodeSetUserData(_cvode_mem, static_cast<void*>(_f.get())));

    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    check_error("CVodeSVtolerances",
                CVodeSVtolerances(_cvode_mem, _reltol, _abstol));

    /* Call CVDense to specify the CVDENSE dense linear solver */
    check_error("CVDense", CVDense(_cvode_mem, _num_equations));

    if (_f->hasJacobian())
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
                    CVDlsSetDenseJacFn(_cvode_mem, df_wrapped));
    }
}

bool CVodeSolverImpl::solve(const double t_end)
{
    realtype t_reached;
    check_error("CVode solve",
                CVode(_cvode_mem, t_end, _y, &t_reached, CV_NORMAL));
    _t = t_reached;

    // check_error asserts that t_end == t_reached and that solving the ODE
    // went fine. Otherwise the program will be aborted. Therefore, we don't
    // have to check manually for errors here and can always savely return true.
    return true;
}

void CVodeSolverImpl::getYDot(const double t, double const* const y,
                              double* const y_dot)
{
    assert(_f != nullptr);
    _f->call(t, y, y_dot);
}

CVodeSolverImpl::~CVodeSolverImpl()
{
    printStats(_cvode_mem);

    N_VDestroy_Serial(_y);
    N_VDestroy_Serial(_abstol);
    CVodeFree(&_cvode_mem);
}

CVodeSolver::CVodeSolver(BaseLib::ConfigTree const& config,
                         unsigned const num_equations)
    : _impl{new CVodeSolverImpl{config, num_equations}}
{
}

void CVodeSolver::setTolerance(const double* abstol, const double reltol)
{
    _impl->setTolerance(abstol, reltol);
}

void CVodeSolver::setTolerance(const double abstol, const double reltol)
{
    _impl->setTolerance(abstol, reltol);
}

void CVodeSolver::setFunction(std::unique_ptr<detail::FunctionHandles>&& f)
{
    _impl->setFunction(std::move(f));
}

void CVodeSolver::setIC(const double t0, double const* const y0)
{
    _impl->setIC(t0, y0);
}

void CVodeSolver::preSolve()
{
    _impl->preSolve();
}

bool CVodeSolver::solve(const double t_end)
{
    return _impl->solve(t_end);
}

double const* CVodeSolver::getSolution() const
{
    return _impl->getSolution();
}

void CVodeSolver::getYDot(const double t, double const* const y,
                          double* const y_dot) const
{
    _impl->getYDot(t, y, y_dot);
}

double CVodeSolver::getTime() const
{
    return _impl->getTime();
}

CVodeSolver::~CVodeSolver() = default;

}  // namespace ODE
}  // namespace MathLib

#endif  // CVODE_FOUND
