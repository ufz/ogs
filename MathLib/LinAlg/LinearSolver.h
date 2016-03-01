/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHLIB_LINEAR_SOLVER_H
#define MATHLIB_LINEAR_SOLVER_H

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{

//! Basic linear solver interface.
template<typename Matrix, typename Vector>
class LinearSolver
{
public:
    /*! Solves the linear equation system \f$ Ax=b \f$ for \f$x \f$
     *
     * \retval true if the system has been successfully solved
     * \retval false otherwise
     */
    virtual bool solve(Matrix& A, Vector& b, Vector& x) = 0;

    virtual ~LinearSolver() = default;
};

// TODO also pass a name argument, or done within config?
/*! Creates a new linear solver instance.
 *
 * \tparam Matrix the matrix type used in the equation system.
 * \tparam Vector the vector type used in the equation system.
 * \tparam Solver the type of the linear solver to be created.
 *
 * \param config configuration options.
 */
template<typename Matrix, typename Vector, typename Solver>
std::unique_ptr<LinearSolver<Matrix, Vector> >
createLinearSolver(BaseLib::ConfigTree const*const config);

/*! Creates a new linear solver instance; this method creates
 *  a linear solver of the default type for the given
 *  \c Matrix and \c Vector types.
 *
 * Said default type is spcified in the LinearSolver.cpp file.
 *
 * This method is an easy way to generate a linear solver, e.g.,
 * in tests.
 *
 * \tparam Matrix the matrix type used in the equation system.
 * \tparam Vector the vector type used in the equation system.
 *
 * \param config configuration options.
 */
template<typename Matrix, typename Vector>
std::unique_ptr<LinearSolver<Matrix, Vector> >
createLinearSolver(BaseLib::ConfigTree const*const config);

} // namespace MathLib

#endif // MATHLIB_LINEAR_SOLVER_H
