/**
 * \file
 *
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
namespace HeatTransportBHE
{
namespace BHE
{
struct RefrigerantParameters
{
    /**
     * dynamics viscosity of the refrigerant
     * unit is kg m-1 sec-1
     */
    double const mu_r;

    /**
     * density of the refrigerant
     * unit is kg m-3
     */
    double const rho_r;

    /**
     * thermal conductivity of the refrigerant
     * unit is kg m sec^-3 K^-1
     */
    double const lambda_r;

    /**
     * specific heat capacity of the refrigerant
     * unit is m^2 sec^-2 K^-1
     */
    double const heat_cap_r;

    /**
     * longitudinal dispersivity of the
     * referigerant flow in the pipeline
     */
    double const alpha_L;
};
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
