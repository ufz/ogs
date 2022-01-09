/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace BaseLib
{
class ConfigTree;
}
namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
struct GroutParameters
{
    /**
     * density of the grout
     * unit is kg m-3
     */
    double const rho_g;

    /**
     * porosity of the grout
     * unit is [-]
     */
    double const porosity_g;

    /**
     * specific heat capacity of the grout
     * unit is m^2 sec^-2 K^-1
     */
    double const heat_cap_g;

    /**
     * thermal conductivity of the grout
     * unit is kg m sec^-3 K^-1
     */
    double const lambda_g;
};

GroutParameters createGroutParameters(BaseLib::ConfigTree const& config);

}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
