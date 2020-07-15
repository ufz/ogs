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

namespace ParameterLib
{
template <typename T>
struct Parameter;
}

namespace ProcessLib::HeatConduction
{
struct HeatConductionProcessData
{
    ParameterLib::Parameter<double> const& thermal_conductivity;
    ParameterLib::Parameter<double> const& heat_capacity;
    ParameterLib::Parameter<double> const& density;

    /// If set mass lumping will be applied to the equation.
    bool const mass_lumping;
};
}  // namespace ProcessLib::HeatConduction
