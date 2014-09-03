/**
 * \author Norihiro Watanabe
 * \date   2012-06-25
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
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

Picard::Picard()
: _normType(VecNormType::INFINITY_N),
  _abs_tol(std::numeric_limits<double>::max()), _rel_tol(1e-6), _max_itr(25),
  _printErrors(false), _n_iterations(0), _abs_error(.0), _rel_error(.0)
{
}

template<class T_FUNCTOR, class T_VALUE>
bool Picard::solve(T_FUNCTOR &functor,  const T_VALUE &x0, T_VALUE &x_new)
{
    T_VALUE x_old(x0);
    T_VALUE dx(x0);

    const bool checkAbsError = (_abs_tol<std::numeric_limits<double>::max());
    const bool checkRelError = (_rel_tol<std::numeric_limits<double>::max());

    INFO("------------------------------------------------------------------");
    INFO("*** PICARD nonlinear solver");
    INFO("-> iteration started");
    bool converged = false;
    std::size_t itr_cnt = 0;
    double x_norm = -1.;
    double abs_error = -1.;
    double rel_error = -1.;
    for (itr_cnt=0; itr_cnt<_max_itr; itr_cnt++) {
        functor(x_old, x_new);
        dx = x_new;
        dx -= x_old;

        abs_error = norm(dx, _normType);
        if (checkRelError) {
            x_norm = norm(x_new, _normType);
            if (x_norm>.0)
                rel_error = abs_error / x_norm;
            else
                rel_error = abs_error;
        }
        converged = (abs_error < _abs_tol && rel_error < _rel_tol);
        if (_printErrors)
            INFO("-> %d: ||dx||=%1.3e, ||x||=%1.3e, ||dx||/||x||=%1.3e", itr_cnt, abs_error, x_norm, rel_error);

#ifdef NDEBUG
        printout(std::cout, itr_cnt, x_new, dx);
#endif
        if (converged) {
            break;
        }
        x_old = x_new;
    }

    INFO("-> iteration finished");
    if (_max_itr==1) {
        INFO("status    : iteration not required");
    } else {
        INFO("status    : %s", (converged ? "CONVERGED" : "***DIVERGED***"));
    }
    INFO("iteration : %d/%d", itr_cnt, _max_itr);
    if (checkAbsError)
        INFO("abs error : %1.3e (tolerance=%1.3e)", abs_error, _abs_tol);
    if (checkRelError)
        INFO("rel error : %1.3e (tolerance=%1.3e)", rel_error, _rel_tol);
    INFO("norm type : %s", convertVecNormTypeToString(_normType).c_str());
    INFO("------------------------------------------------------------------");

    this->_n_iterations = itr_cnt;
    this->_abs_error = abs_error;
    this->_rel_error = rel_error;

    return converged;
}


#ifdef NDEBUG
template<class T_VALUE>
inline void Picard::printout(std::ostream& os, std::size_t i, T_VALUE& x_new, T_VALUE& dx)
{
    os << "-> " << i <<": x=(";
    for (std::size_t i=0; i<x_new.size(); i++)
        os << x_new[i] << " ";
    os << "), dx=(";
    for (std::size_t i=0; i<dx.size(); i++)
        os << dx[i] << " ";
    os << ")\n";
}

// in case of double
template<>
inline void Picard::printout(std::ostream& os, std::size_t i, double& x_new, double& dx)
{
    os << "-> " << i <<": x=" << x_new << ", dx=" << dx << "\n";
}
#endif

}

} //end
