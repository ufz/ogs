/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TESAssemblyParams.h"

namespace ProcessLib
{
namespace TES
{
class TESFEMReactionAdaptor;

struct TESLocalAssemblerData
{
    TESLocalAssemblerData(AssemblyParams const& ap_, const unsigned num_int_pts,
                          const unsigned dimension);

    ~TESLocalAssemblerData();

    AssemblyParams const& ap;

    // integration point quantities
    std::vector<double> solid_density;
    std::vector<double> reaction_rate;  // dC/dt * rho_SR_dry
    std::vector<std::vector<double>>
        velocity;  // vector of velocities for each integration point

    // integration point values of unknowns -- temporary storage
    double p = std::numeric_limits<double>::quiet_NaN();  // gas pressure
    double T = std::numeric_limits<double>::quiet_NaN();  // temperature
    double vapour_mass_fraction = std::numeric_limits<
        double>::quiet_NaN();  // fluid mass fraction of the second component

    // temporary storage for some properties
    // values do not change during the assembly of one integration point
    double rho_GR = std::numeric_limits<double>::quiet_NaN();
    double p_V =
        std::numeric_limits<double>::quiet_NaN();  // vapour partial pressure
    double qR = std::numeric_limits<
        double>::quiet_NaN();  // reaction rate, use this in assembly!!!

    std::unique_ptr<TESFEMReactionAdaptor> const reaction_adaptor;

    // variables at previous timestep
    std::vector<double> solid_density_prev_ts;
    std::vector<double> reaction_rate_prev_ts;  // could also be calculated from
                                                // solid_density_prev_ts
};
}
}  // namespaces
