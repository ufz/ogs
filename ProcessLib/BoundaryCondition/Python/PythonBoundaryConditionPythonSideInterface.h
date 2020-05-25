/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
//! Base class for boundary conditions.
//! This class will get Python bindings and is intended to be to be derived in
//! Python.
class PythonBoundaryConditionPythonSideInterface
{
public:
    /*!
     * Computes Dirichlet boundary condition values for the provided arguments
     * (time, position of the node, node id, primary variables at the node).
     *
     * \return a pair (is_dirichlet, value) indicating if a Dirichlet BC shall
     * be set at that node and (if so) the value of the Dirichlet BC at that
     * node.
     */
    virtual std::pair<bool, double> getDirichletBCValue(
        double /*t*/, std::array<double, 3> /*x*/, std::size_t /*node_id*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        overridden_essential_ = false;
        return {false, std::numeric_limits<double>::quiet_NaN()};
    }

    /*!
     * Computes the flux for the provided arguments (time, position, primary
     * variables at that position).
     *
     * \return a pair (is_natural, flux, flux_jacobian) indicating if a natural
     * BC shall be set at that position and (if so) the flux at that node and
     * the derivative of the flux w.r.t. all primary variables.
     */
    virtual std::tuple<bool, double, std::vector<double>> getFlux(
        double /*t*/,
        std::array<double, 3> /*x*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        overridden_natural_ = false;
        return std::tuple<bool, double, std::vector<double>>{
            false, std::numeric_limits<double>::quiet_NaN(), {}};
    }

    //! Tells if getDirichletBCValue() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre getDirichletBCValue() must already have been called
    //! once.
    bool isOverriddenEssential() const { return overridden_essential_; }

    //! Tells if getFlux() has been overridden in the derived class in Python.
    //!
    //! \pre getFlux() must already have been called once.
    bool isOverriddenNatural() const { return overridden_natural_; }

    virtual ~PythonBoundaryConditionPythonSideInterface() = default;

private:
    //! Tells if getDirichletBCValue() has been overridden in the derived class
    //! in Python.
    mutable bool overridden_essential_ = true;
    //! Tells if getFlux() has been overridden in the derived class in Python.
    mutable bool overridden_natural_ = true;
};
}  // namespace ProcessLib
