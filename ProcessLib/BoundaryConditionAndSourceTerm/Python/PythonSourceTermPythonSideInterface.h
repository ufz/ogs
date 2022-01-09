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
    virtual std::pair<double, std::vector<double>> getFlux(
        double /*t*/, std::array<double, 3> const& /*x*/,
        std::vector<double> const& /*primary_variables*/) const = 0;

    virtual ~PythonSourceTermPythonSideInterface() = default;
};
}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
