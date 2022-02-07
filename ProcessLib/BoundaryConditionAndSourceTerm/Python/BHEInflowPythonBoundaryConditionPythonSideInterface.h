/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
        _overridden_essential = false;
        return std::tuple<double,
                          std::vector<double>,
                          std::vector<double>,
                          std::vector<int>,
                          std::vector<double>>{
            std::numeric_limits<double>::quiet_NaN(), {}, {}, {}, {}};
    }

    /*!
     * transfer BHE network dataframe to TESPy and get Tin and flow rate from
     * TESPy
     *
     * \return a tuple (if use tespyThermalSolver, if convergence achieved
     * in tespy, BHE Tin value, BHE flow rate from TESPy)
     * indicating if tespyThermalSolver shall be used at that position, if
     * themal convergence has been achieved in the tespy, the new
     * inflow temperature and flow rate for all BHEs.
     */
    virtual std::tuple<bool, bool, std::vector<double>, std::vector<double>>
    tespySolver(double /*t*/,
                std::vector<double> const& /*Tin_val*/,
                std::vector<double> const& /*Tout_val*/) const
    {
        _overridden_tespy = false;
        return std::tuple<bool, bool, std::vector<double>, std::vector<double>>{
            false, false, {}, {}};
    }

    /*!
     * Transfer BHE network dataframe to python script where function
     * serverCommunication takes t, dt, Tin_val, Tout_val, BHE_flowrate as
     * arguments and returns updated Tin_val and flowrate which will be used
     * as (python) boundary conditions for the next time step.
     * The communication happens in preTimestep so that the updated values
     * are used in the current time step.
     * Specific use case: Interface between OGS and SimulationX via TCP/IP
     * where OGS is connected as one client
     *
     * \return a tuple (BHE Tin value, BHE flow rate for all BHEs from
     * serverCommunication)
     */
    virtual std::tuple<std::vector<double>, std::vector<double>>
    serverCommunication(double /*t*/,
                        double /*dt*/,
                        std::vector<double> const& /*Tin_val*/,
                        std::vector<double> const& /*Tout_val*/,
                        std::vector<double> const& /*BHE_flowrate*/) const
    {
        _overridden_server_communication = false;
        return {};
    }

    //! Tells if initializeDataContainer() has been overridden in the derived
    //! class in Python.
    //!
    //! \pre initializeDataContainer() must already have been called
    //! once.
    bool isOverriddenEssential() const { return _overridden_essential; }

    //! Tells if tespySolver() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre tespySolver() must already have been called once.
    bool isOverriddenTespy() const { return _overridden_tespy; }

    //! Tells if serverCommunication() has been overridden in the derived class
    //! in Python.
    //!
    //! \pre serverCommunication() must already have been called once.
    bool isOverriddenServerCommunication() const
    {
        return _overridden_server_communication;
    }

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
    mutable bool _overridden_essential = true;
    //! Tells if tespySolver() has been overridden in the derived class
    //! in Python.
    mutable bool _overridden_tespy = true;
    //! Tells if serverCommunication() has been overridden in the derived class
    //! in Python.
    mutable bool _overridden_server_communication = true;
};
}  // namespace ProcessLib
