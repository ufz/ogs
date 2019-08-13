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

template <typename T>
struct Parameter;

namespace GroundwaterFlow
{
struct GroundwaterFlowProcessData final
{
    ParameterLib::Parameter<double> const& hydraulic_conductivity;
};

} // namespace GroundwaterFlow
} // namespace ProcessLib
