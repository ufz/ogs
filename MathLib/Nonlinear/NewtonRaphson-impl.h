/**
 * \author Norihiro Watanabe
 * \date   2012-06-25
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <limits>

#include "logog/include/logog.hpp"

namespace MathLib
{

namespace Nonlinear
{

NewtonRaphson::NewtonRaphson()
: _normType(VecNormType::INFINITY_N),
  _r_abs_tol(std::numeric_limits<double>::max()), _r_rel_tol(1e-6), _dx_rel_tol(.0),
  _max_itr(25), _printErrors(true), _n_iterations(0), _r_abs_error(.0), _r_rel_error(.0), _dx_rel_error(.0)
{
}

template<class F_RESIDUAL, class F_DX, class T_VALUE>
bool NewtonRaphson::solve(F_RESIDUAL &f_residual, F_DX &f_dx, const T_VALUE &x0, T_VALUE &x_new)
{
    T_VALUE r(x0), dx(x0);
    r = .0;
    dx = .0;

    const bool checkAbsResidualError = (_r_abs_tol<std::numeric_limits<double>::max());
    const bool checkRelResidualError = (_r_rel_tol<std::numeric_limits<double>::max());
    const bool checkRelDxError = (_dx_rel_tol>.0);
    const bool needXNorm =  (checkRelResidualError || checkRelDxError);

    INFO("------------------------------------------------------------------");
    INFO("*** NEWTON-RAPHSON nonlinear solver");
    INFO("-> iteration started");

    double x_norm = -1.;
    double dx_norm = -1.;
    double r_norm = -1.;
    bool converged = false;

    // evaluate initial residual
    x_new = x0;
    f_residual(x_new, r);
    // check convergence
    r_norm = norm(r, _normType);
    dx_norm = norm(dx, _normType);
    if (needXNorm)
        x_norm = norm(x_new, _normType);
    converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
    if (_printErrors)
        INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", 0, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

    std::size_t itr_cnt = 0;
    if (!converged) {
        for (itr_cnt=1; itr_cnt<_max_itr; itr_cnt++) {
            // solve dx=-J^-1*r
            f_dx(x_new, r, dx);
            x_new += dx;
            // evaluate residual
            f_residual(x_new, r);
#ifdef NDEBUG
            printout(std::cout, itr_cnt, x_new, r, dx);
#endif
            // check convergence
            r_norm = norm(r, _normType);
            dx_norm = norm(dx, _normType);
            if (needXNorm)
                x_norm = norm(x_new, _normType);
            converged = ((r_norm < _r_abs_tol && r_norm < _r_rel_tol*x_norm)
                        || (checkRelDxError && dx_norm < _dx_rel_tol*x_norm));
            if (_printErrors)
                INFO("-> %d: ||r||=%1.3e, ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, r_norm, dx_norm, x_norm, x_norm==0 ? dx_norm : dx_norm/x_norm);

            if (converged)
                break;
        }
    }

    INFO("-> iteration finished");
    if (_max_itr==1) {
        INFO("status    : iteration not required");
    } else {
        INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
    }
    INFO("iteration : %d/%d", itr_cnt, _max_itr);
    if (checkAbsResidualError)
        INFO("abs. res. : %1.3e (tolerance=%1.3e)", r_norm, _r_abs_tol);
    if (checkRelResidualError)
        INFO("rel. res. : %1.3e (tolerance=%1.3e)", x_norm==0?r_norm:r_norm/x_norm, _r_rel_tol);
    if (checkRelDxError)
        INFO("dx        : %1.3e (tolerance=%1.3e)", x_norm==0?dx_norm:dx_norm/x_norm, _dx_rel_tol);
    INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
    INFO("------------------------------------------------------------------");

    this->_n_iterations = itr_cnt;
    this->_r_abs_error = r_norm;
    this->_r_rel_error = (x_norm==0 ? r_norm : r_norm / x_norm);
    this->_dx_rel_error = (x_norm==0 ? dx_norm : dx_norm / x_norm);

    return converged;
}


#ifdef NDEBUG
template<class T_VALUE>
inline void NewtonRaphson::printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& r, T_VALUE& dx)
{
    os << "-> " << i <<": x=(";
    for (std::size_t i=0; i<x_new.size(); i++)
        os << x_new[i] << " ";
    os << "), r=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << r[i] << " ";
    os << "), dx=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << dx[i] << " ";
    os << ")\n";
}

// in case of double
template<>
inline void NewtonRaphson::printout(std::ostream& os, std::size_t i, double& x_new, double& r, double& dx)
{
    os << "-> " << i <<": x=" << x_new << ", r=" << r << ", dx=" << dx << "\n";
}
#endif

}

} //end
