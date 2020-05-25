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
//! Base class for BHENetwork.
//! This class will get Python bindings and is intended to be to be derived in
//! Python.
class BHEInflowPythonBoundaryConditionPythonSideInterface
{
public:
    /*!
     * Initialize network dataframe
     * return a tuple (time, BHE inflow temperature, BHE outflow
     * temperature, BHE outflow bc node id, BHE flowrate)
     * set at that position and the parameters of the BHE network.
     */
    virtual std::tuple<double /*time*/,
                       std::vector<double> /*Tin_val*/,
                       std::vector<double> /*Tout_val*/,
                       std::vector<int> /*bc_out_ids*/,
                       std::vector<double> /*BHE_flowrate*/>
    initializeDataContainer() const
    {
        overridden_essential_ = false;
        return std::tuple<double,
                          std::vector<double>,
                          std::vector<double>,
                          std::vector<int>,
                          std::vector<double>>{
            std::numeric_limits<double>::quiet_NaN(), {}, {}, {}, {}};
    }

    /*!
     * transfer BHE network dataframe to TESPy and get Tin from TESPy
     *
     * \return a tuple (if use tespyThermalSolver, if convergence achieved
     * in tespy, BHE Tin and Tout value from TESPy)
     * indicating if tespyThermalSolver shall be used at that position, if
     * themal convergence has been achieved in the tespy and the new
     * inflow temperature of all BHEs.
     */
    virtual std::tuple<bool, bool, std::vector<double>> tespyThermalSolver(
        double /*t*/,
        std::vector<double> const& /*Tin_val*/,
        std::vector<double> const& /*Tout_val*/) const
    {
        overridden_tespyThermal_ = false;
        return std::tuple<bool, bool, std::vector<double>>{false, false, {}};
    }

    /*!
     * call Tespy hydraulic solver to get flow velocity in each pipe
     *
     * \return a tuple (if use tespyHydroSolver,  f_velocity) indicating if
     * tespyHydroSolver shall be used at that position and the flow velocity
     * in each pipe of all BHEs.
     */
    virtual std::tuple<bool, std::vector<double>> tespyHydroSolver(
        double /*t*/) const
    {
        overridden_tespyHydro_ = false;
        return std::tuple<bool, std::vector<double>>{false, {}};
    }

    //! Tells if initializeDataContainer() has been overridden in the derived
    //! class in Python.
    //!
    //! \pre initializeDataContainer() must already have been called
    //! once.
    bool isOverriddenEssential() const { return overridden_essential_; }

    //! Tells if tespyThermalSolver() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre tespyThermalSolver() must already have been called once.
    bool isOverriddenTespyThermal() const { return overridden_tespyThermal_; }

    //! Tells if tespyHydroSolver() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre tespyHydroSolver() must already have been called once.
    bool isOverriddenTespyHydro() const { return overridden_tespyHydro_; }

    // BHE network dataframe container
    std::tuple<double,
               std::vector<double>,
               std::vector<double>,
               std::vector<int>,
               std::vector<double>>
        dataframe_network;

    virtual ~BHEInflowPythonBoundaryConditionPythonSideInterface() = default;

private:
    //! Tells if initializeDataContainer() has been overridden in the derived
    //! class in Python.
    mutable bool overridden_essential_ = true;
    //! Tells if tespyThermalSolver() has been overridden in the derived class
    //! in Python.
    mutable bool overridden_tespyThermal_ = true;
    //! Tells if tespyHydroSolver() has been overridden in the derived class in
    //! Python.
    mutable bool overridden_tespyHydro_ = true;
};
}  // namespace ProcessLib
