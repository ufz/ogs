/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
//! Base class for BHENetwork.
//! This class will get Python bindings and is intended to be to be derived in
//! Python.
class BHEInflowPythonBoundaryConditionPythonSideInterface
{
public:
    /*!
     * Initialize network dataframe
     * return a tuple (is_natural, BHE inflow temperature, BHE outflow
     * temperature, BHE outflow bc node id, time) indicating if a natural BC shall be
     * set at that position and the parameters of the BHE network.
     */
    virtual std::tuple<bool,
                       std::vector<double> /*Tin_val*/,
                       std::vector<double> /*Tout_val*/,
                       std::vector<int> /*bc_out_ids*/,
                       double /*time*/>
        initializeDataContainer() const
    {
        _overridden_essential = false;
        return std::tuple<bool,
                          std::vector<double>,
                          std::vector<double>,
                          std::vector<int>,
                          double>{
            false, {}, {}, {}, std::numeric_limits<double>::quiet_NaN()};
    }

    /*!
     * transfer BHE network dataframe to TESPy and get Tin from TESPy
     *
     * \return a tuple (is_natural, if convergence, Tin_val from TESPy) indicating if a
     * natural BC shall be set at that position and the new inflow temperature
     * of all BHEs
     */
    virtual std::tuple<bool, bool, std::vector<double>> tespyThermalSolver(
        double /*t*/,
        std::vector<double> const& /*Tin_val*/,
        std::vector<double> const& /*Tout_val*/) const
    {
        _overridden_natural = false;
        return std::tuple<bool, bool, std::vector<double>>{false, false, {}};
    }

    /*!
     * call Tespy hydraulic solver to get flow velocity in each pipe
     *
     * \return a tuple (is_natural,  f_velocity ) indicating if a
     * natural BC shall be set at that position and the flow velocity in each
     * pipe of all BHEs
     */
    virtual std::tuple<bool, std::vector<double>> tespyHydroSolver() const
    {
        _overridden_natural = false;
        return std::tuple<bool, std::vector<double>>{false, {}};
    }

    //! Tells if initializeDataContainer() has been overridden in the derived
    //! class in Python.
    //!
    //! \pre initializeDataContainer() must already have been called
    //! once.
    bool isOverriddenEssential() const { return _overridden_essential; }

    //! Tells if tespySolver() has been overridden in the derived class in
    //! Python.
    //!
    //! \pre tespySolver() must already have been called once.
    bool isOverriddenNatural() const { return _overridden_natural; }

    // BHE network dataframe container
    std::tuple<bool, std::vector<double>, std::vector<double>, std::vector<int>, double>
        dataframe_network;

    virtual ~BHEInflowPythonBoundaryConditionPythonSideInterface() = default;

private:
    //! Tells if initializeDataContainer() has been overridden in the derived
    //! class in Python.
    mutable bool _overridden_essential = true;
    //! Tells if tespySolver() has been overridden in the derived class in
    //! Python.
    mutable bool _overridden_natural = true;
};
}  // namespace ProcessLib
