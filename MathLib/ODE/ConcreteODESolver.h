/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "FunctionHandles.h"
#include "ODESolver.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
namespace ODE
{
template <unsigned NumEquations>
std::unique_ptr<ODESolver<NumEquations>> createODESolver(
    BaseLib::ConfigTree const& config);

//! \addtogroup ExternalODESolverInterface
//! @{

/*! ODE solver with a bounds-safe interface using a concrete implementation.
 *
 * This class makes contact between the abstract \c ODESolver interface and a
 * certain solver \c Implementation.
 *
 * The interface of this class inherits the array bounds checking from \c
 * ODESolver.
 * Its methods forward calls to the \c Implementation erasing array bounds info
 * by passing \c std::array and Eigen::Matrix arguments as raw pointers.
 *
 * Thus the \c Implementation does not need template to be templated
 * which makes it possible for \c Implementation to use the pimpl idiom whereby
 * the headers of external libraries only have to be included in the final cpp
 * file. Thereby our namespaces do not get polluted with symbols from external
 * libraries.
 *
 * \tparam NumEquations the number of equations in the ODE system.
 * \tparam Implementation a concrete ODE solver implementation used as a
 * backend.
 */
template <typename Implementation, unsigned NumEquations>
class ConcreteODESolver final : public ODESolver<NumEquations>,
                                private Implementation
{
public:
    void setFunction(Function<NumEquations> f,
                     JacobianFunction<NumEquations> df) override
    {
        Implementation::setFunction(
            std::unique_ptr<detail::FunctionHandlesImpl<NumEquations>>{
                new detail::FunctionHandlesImpl<NumEquations>{f, df}});
    }

    void setTolerance(const std::array<double, NumEquations>& abstol,
                      const double reltol) override
    {
        Implementation::setTolerance(abstol.data(), reltol);
    }

    void setTolerance(const double abstol, const double reltol) override
    {
        Implementation::setTolerance(abstol, reltol);
    }

    void setIC(const double t0,
               std::initializer_list<double> const& y0) override
    {
        assert(y0.size() == NumEquations);
        Implementation::setIC(t0, y0.begin());
    }

    void setIC(const double t0,
               Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> const&
                   y0) override
    {
        Implementation::setIC(t0, y0.data());
    }

    void preSolve() override { Implementation::preSolve(); }
    bool solve(const double t) override { return Implementation::solve(t); }
    MappedConstVector<NumEquations> getSolution() const override
    {
        return MappedConstVector<NumEquations>{Implementation::getSolution()};
    }
    double getTime() const override { return Implementation::getTime(); }
    Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> getYDot(
        const double t, const MappedConstVector<NumEquations>& y) const override
    {
        Eigen::Matrix<double, NumEquations, 1, Eigen::ColMajor> y_dot;
        Implementation::getYDot(t, y.data(), y_dot.data());
        return y_dot;
    }

private:
    //! Instances of this class shall only be constructed by createODESolver().
    ConcreteODESolver(BaseLib::ConfigTree const& config)
        : Implementation{config, NumEquations}
    {
    }

    friend std::unique_ptr<ODESolver<NumEquations>>
    createODESolver<NumEquations>(BaseLib::ConfigTree const& config);
};

//! @}

}  // namespace ODE
}  // namespace MathLib
