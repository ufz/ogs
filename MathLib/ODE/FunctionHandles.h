/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ODESolverTypes.h"

namespace MathLib
{
namespace ODE
{
namespace detail
{
//! \addtogroup ExternalODESolverInterface
//! @{

/*! Interface providing acces to functions computing \f$\dot y\f$ and its
 *  Jacobian to code interfacing with external ODE solver libraries.
 *
 * \note
 * The methods of this class accept raw pointers as arguments, not Eigen::Matrix
 * types.
 */
class FunctionHandles
{
public:
    //! Calls a function computing \f$\dot y\f$.
    //! \returns true or false indicating whether the function succeeded.
    virtual bool call(const double t, double const* const y,
                      double* const ydot) = 0;

    //! Calls a function computing \f$\mathtt{jac} := \partial \dot y/\partial
    //! y\f$.
    //! \returns true or false indicating whether the function succeeded.
    virtual bool callJacobian(const double t,
                              double const* const y,
                              double* const ydot,
                              double* const jac) = 0;

    //! Tells whether a Jacobian function has been set.
    virtual bool hasJacobian() const = 0;

    //! Returns the number of equations in the ODE system.
    virtual unsigned getNumberOfEquations() const = 0;

    virtual ~FunctionHandles() = default;
};

//! Function handles for an ODE system of \c N equations.
template <unsigned N>
struct FunctionHandlesImpl final : public FunctionHandles
{
    FunctionHandlesImpl(Function<N>& f, JacobianFunction<N>& df) : f(f), df(df)
    {
    }

    /*! Calls the stored function \c f computing \f$\dot y\f$.
     *
     * The raw pointers passed to this method are wrapped in some Eigen::Map
     * objects before being passed to \c f. Thereby the information about the
     * size of the vectors is restored. No memory is copied for that.
     *
     * \returns true or false indicating whether the function succeeded.
     */
    bool call(const double t, const double* const y,
              double* const ydot) override
    {
        if (f)
        {
            MappedVector<N> ydot_mapped{ydot};
            return f(t, MappedConstVector<N>{y}, ydot_mapped);
        }
        return false;
    }

    /*! Calls the stored function computing
     *  \f$\mathtt{jac} := \partial \dot y/\partial y\f$.
     *
     * \returns true or false indicating whether the function succeeded.
     * \see call()
     */
    bool callJacobian(const double t, const double* const y, double* const ydot,
                      double* const jac) override
    {
        if (df)
        {
            MappedMatrix<N, N> jac_mapped{jac};
            return df(t,
                      MappedConstVector<N>{y},
                      MappedConstVector<N>{ydot},
                      jac_mapped);
        }
        return false;
    }

    bool hasJacobian() const override { return df != nullptr; }
    unsigned getNumberOfEquations() const override { return N; }
    Function<N> f;
    JacobianFunction<N> df;
};

//! @}

}  // namespace detail
}  // namespace ODE
}  // namespace MathLib
