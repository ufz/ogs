/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <Eigen/Core>

namespace MathLib
{
namespace ODE
{
//! \addtogroup ExternalODESolverInterface
//! @{

/*! This type can be used like an \f$N \times M\f$ Eigen::Matrix, but it
 *  does not manage the storage for its entries on its own.
 *
 * \tparam N number of rows
 * \tparam M number of columns
 *
 * \see https://eigen.tuxfamily.org/dox/classEigen_1_1Map.html
 */
template <int N, int M>
using MappedMatrix = Eigen::Map<Eigen::Matrix<double, N, M, Eigen::ColMajor>>;

//! Behaves like a \c const Eigen::Matrix. \see MappedMatrix
template <int N, int M>
using MappedConstMatrix =
    Eigen::Map<const Eigen::Matrix<double, N, M, Eigen::ColMajor>>;

//! \see MappedMatrix
template <int N>
using MappedVector = MappedMatrix<N, 1>;

//! \see MappedConstMatrix
template <int N>
using MappedConstVector = MappedConstMatrix<N, 1>;

/*! A function computing \c ydot as a function of \c t and \c y.
 *
 *  The function returns true or false indecating whether it succeeded.
 *
 * \tparam N the number of equations in the ODE system.
 */
template <unsigned N>
using Function = std::function<bool(
    const double t, MappedConstVector<N> const& y, MappedVector<N>& ydot)>;

/*! A function computing \f$\mathtt{jac} := \partial \dot y/\partial y\f$.
 *
 * The function returns true or false indecating whether it succeeded.
 *
 * \tparam N the number of equations in the ODE system.
 */
template <unsigned N>
using JacobianFunction = std::function<bool(const double t,
                                            MappedConstVector<N> const& y,
                                            MappedConstVector<N> const& ydot,
                                            MappedMatrix<N, N>& jac)>;

//! @}

}  // namespace ODE
}  // namespace MathLib
