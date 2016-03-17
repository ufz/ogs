/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

extern "C" {
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
}

#include "CVodeSolver.h"

#include <cassert>
#include <logog/include/logog.hpp>

namespace
{
void printStats(void* cvode_mem)
{
	long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
	int flag;

	flag = CVodeGetNumSteps(cvode_mem, &nst);
	// check_flag(&flag, "CVodeGetNumSteps", 1);
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	// check_flag(&flag, "CVodeGetNumRhsEvals", 1);
	flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
	// check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
	flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
	// check_flag(&flag, "CVodeGetNumErrTestFails", 1);
	flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
	// check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
	flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
	// check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

	flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
	// check_flag(&flag, "CVDlsGetNumJacEvals", 1);
	flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
	// check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

	flag = CVodeGetNumGEvals(cvode_mem, &nge);
	// check_flag(&flag, "CVodeGetNumGEvals", 1);

	DBUG("\nFinal Statistics:");
	DBUG("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld",
	     nst, nfe, nsetups, nfeLS, nje);
	DBUG("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n", nni, ncfn, netf,
	     nge);
}
}

namespace MathLib
{
/**
 * This class provides concrete access to Sundials' CVode solver.
 *
 * OdeSolver (implicitly bounds checked, agnostic to concrete ODE solver)
 *  |
 *  | Dynamic polymorphism
 *  v
 * ConcreteOdeSolver (implicitly bounds checked, interfaces with a specific
 * library)
 *  |
 *  | Forward calls, disable bounds checking, no need for templates anymore
 *  v
 * Implementation = CVodeSolverInternal (no templates)
 *  |
 *  | Pimpl (hide implementation, do not include 3rd party libs in header)
 *  v
 * CVodeSolverImpl
 *
 * This scheme might be a general way for accessing 3rd party libraries.
 */
class CVodeSolverImpl
{
	static_assert(std::is_same<realtype, double>::value,
	              "cvode's realtype is not the same as double");

public:
	CVodeSolverImpl(BaseLib::ConfigTree const& config);
	void init(const unsigned num_equations);

	void setTolerance(const double* abstol, const double reltol);
	void setTolerance(const double abstol, const double reltol);

	void setFunction(std::unique_ptr<FunctionHandles>&& f);

	void setIC(const double t0, double const* const y0);

	void preSolve();
	void solve(const double t_end);

	double const* getSolution() const { return NV_DATA_S(_y); }
	double getTime() const { return _t; }
	bool getYDot(const double t, double const* const y, double* const ydot);

	~CVodeSolverImpl();

private:
	N_Vector _y = nullptr;

	realtype _t;

	N_Vector _abstol = nullptr;
	realtype _reltol;

	unsigned _num_equations;
	void* _cvode_mem;

	std::unique_ptr<FunctionHandles> _f;

	int _linear_multistep_method = CV_ADAMS;
	int _nonlinear_solver_iteration = CV_FUNCTIONAL;
};

CVodeSolverImpl::CVodeSolverImpl(const BaseLib::ConfigTree& config)
{
	if (auto const param =
	        config.getConfParamOptional<std::string>("linear_multistep_method"))
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
			ERR("unknown linear multistep method: %s", param->c_str());
			std::abort();
		}
	}

	if (auto const param = config.getConfParamOptional<std::string>(
	        "nonlinear_solver_iteration"))
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
			ERR("unknown nonlinear solver iteration: %s", param->c_str());
			std::abort();
		}
	}
}

void CVodeSolverImpl::init(const unsigned num_equations)
{
	_y = N_VNew_Serial(num_equations);
	_abstol = N_VNew_Serial(num_equations);
	_num_equations = num_equations;

	_cvode_mem =
	    CVodeCreate(_linear_multistep_method, _nonlinear_solver_iteration);

	assert(_cvode_mem != nullptr && _y != nullptr && _abstol != nullptr);
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

void CVodeSolverImpl::setFunction(std::unique_ptr<FunctionHandles>&& f)
{
	_f = std::move(f);
	assert(_num_equations == _f->getNumEquations());

	auto f_wrapped = [](const realtype t, const N_Vector y, N_Vector ydot,
	                    void* function_handles) -> int
	{
		bool successful = static_cast<FunctionHandles*>(function_handles)
		                      ->call(t, NV_DATA_S(y), NV_DATA_S(ydot));
		return successful ? 0 : 1;
	};

	// TODO: check not run twice! move this call somewhere else
	int flag = CVodeInit(_cvode_mem, f_wrapped, 0.0, _y);
	(void)flag;
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

	// int flag = CVodeInit(_cvode_mem, f_wrapped, _t, _y); // TODO: consider
	// CVodeReInit()!
	int flag =
	    CVodeReInit(_cvode_mem, _t, _y);  // TODO: consider CVodeReInit()!
	// if (check_flag(&flag, "CVodeInit", 1)) return(1);

	flag = CVodeSetUserData(_cvode_mem, static_cast<void*>(_f.get()));

	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(_cvode_mem, _reltol, _abstol);
	// if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	/* Call CVDense to specify the CVDENSE dense linear solver */
	flag = CVDense(_cvode_mem, _num_equations);
	// if (check_flag(&flag, "CVDense", 1)) return(1);

	if (_f->hasJacobian())
	{
		auto df_wrapped = [](
		    const long /*N*/, const realtype t, const N_Vector y,
		    const N_Vector ydot, const DlsMat jac, void* function_handles,
		    N_Vector /*tmp1*/, N_Vector /*tmp2*/, N_Vector /*tmp3*/
		    ) -> int
		{
			// Caution: by calling the DENSE_COL() macro we assume that matrices
			//          are stored contiguously in memory!
			bool successful =
			    static_cast<FunctionHandles*>(function_handles)
			        ->callJacobian(t, NV_DATA_S(y), NV_DATA_S(ydot),
			                       DENSE_COL(jac, 0));
			return successful ? 0 : 1;
		};

		flag = CVDlsSetDenseJacFn(_cvode_mem, df_wrapped);
		// if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);
	}
}

void CVodeSolverImpl::solve(const double t_end)
{
	realtype t_reached;
	int flag = CVode(_cvode_mem, t_end, _y, &t_reached, CV_NORMAL);
	_t = t_reached;
	// std::cout << "result at time " << t << " is " << NV_Ith_S(y,0) <<
	// std::endl;
	if (flag != CV_SUCCESS)
	{
		// std::cerr << "ERROR at " << __FUNCTION__ << ":" << __LINE__ <<
		// std::endl;
	}
}

bool CVodeSolverImpl::getYDot(const double t, double const* const y,
                              double* const ydot)
{
	if (_f != nullptr)
	{
		return _f->call(t, y, ydot);
	}

	return false;
}

CVodeSolverImpl::~CVodeSolverImpl()
{
	printStats(_cvode_mem);

	if (_y)
	{
		N_VDestroy_Serial(_y);
		N_VDestroy_Serial(_abstol);
	}

	if (_cvode_mem)
	{
		CVodeFree(&_cvode_mem);
	}
}

CVodeSolverInternal::CVodeSolverInternal(BaseLib::ConfigTree const& config)
    : _impl{new CVodeSolverImpl{config}}
{
}

void CVodeSolverInternal::init(const unsigned num_equations)
{
	_impl->init(num_equations);
}

void CVodeSolverInternal::setTolerance(const double* abstol,
                                       const double reltol)
{
	_impl->setTolerance(abstol, reltol);
}

void CVodeSolverInternal::setTolerance(const double abstol, const double reltol)
{
	_impl->setTolerance(abstol, reltol);
}

void CVodeSolverInternal::setFunction(std::unique_ptr<FunctionHandles>&& f)
{
	_impl->setFunction(std::move(f));
}

void CVodeSolverInternal::setIC(const double t0, double const* const y0)
{
	_impl->setIC(t0, y0);
}

void CVodeSolverInternal::preSolve()
{
	_impl->preSolve();
}

void CVodeSolverInternal::solve(const double t_end)
{
	_impl->solve(t_end);
}

double const* CVodeSolverInternal::getSolution() const
{
	return _impl->getSolution();
}

bool CVodeSolverInternal::getYDot(const double t, double const* const y,
                                  double* const ydot) const
{
	return _impl->getYDot(t, y, ydot);
}

double CVodeSolverInternal::getTime() const
{
	return _impl->getTime();
}

CVodeSolverInternal::~CVodeSolverInternal() = default;

}  // namespace MathLib
