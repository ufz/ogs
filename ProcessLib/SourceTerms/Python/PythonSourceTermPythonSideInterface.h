/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
//! Base class for source terms.
//! This class will get Python bindings and is intended to be to be derived in
//! Python.
class PythonSourceTermPythonSideInterface
{
public:
    /*!
     * Computes the flux for the provided arguments (time, position of the node,
     * primary variables at the node).
     *
     * \return flux Flux of the source term at that node and derivative of the
     * flux w.r.t. all primary variables.
     */
    virtual std::pair<double, std::array<double, 3>> getFlux(
        double /*t*/, std::array<double, 3> /*x*/,
        std::vector<double> const& /*primary_variables*/) const
    {
        return {std::numeric_limits<double>::quiet_NaN(), {}};
    }

    //! Tells if getFlux() has been overridden in the derived class in Python.
    //!
    //! \pre getFlux() must already have been called once.
    bool isOverriddenGetFlux() const { return _overridden_get_flux; }

    virtual ~PythonSourceTermPythonSideInterface() = default;

private:
    //! Tells if getFlux() has been overridden in the derived class in Python.
    mutable bool _overridden_get_flux = true;
};
}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
