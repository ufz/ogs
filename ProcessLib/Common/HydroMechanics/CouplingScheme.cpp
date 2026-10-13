// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CouplingScheme.h"

#include <cmath>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "BaseLib/Logging.h"

namespace ProcessLib::Common::HydroMechanics
{
CouplingScheme parseCouplingScheme(BaseLib::ConfigTree const& config)
{
    auto const coupling_scheme_config =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme}
        config.getConfigSubtreeOptional("coupling_scheme");

    if (!coupling_scheme_config)
    {
        return Monolithic{};
    }

    auto const coupling_scheme_type =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__type}
        coupling_scheme_config->getConfigParameterOptional<std::string>("type");

    if (!coupling_scheme_type)
    {
        OGS_FATAL(
            "The coupling scheme configuration is missing its <type> "
            "element.\n"
            "Use the nested form:\n"
            "    <coupling_scheme>\n"
            "        <type>staggered</type>\n"
            "    </coupling_scheme>\n"
            "with the type named 'monolithic' or 'staggered'.");
    }

    auto const& type = *coupling_scheme_type;

    if (type == "monolithic")
    {
        return Monolithic{};
    }

    if (type != "staggered")
    {
        OGS_FATAL(
            "Unknown coupling scheme type '{:s}'. Valid options are "
            "'monolithic' and 'staggered'.",
            type);
    }

    // Default value is 0.5, which is recommended by [Mikelic & Wheeler].
    double const fixed_stress_stabilization_parameter =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_stabilization_parameter}
        coupling_scheme_config->getConfigParameter<double>(
            "fixed_stress_stabilization_parameter", 0.5);

    DBUG("Using value {:g} for coupling parameter of staggered scheme.",
         fixed_stress_stabilization_parameter);

    {  // Check parameter value. Optimum is not a-priori known, but within
       // certain interval [Storvik & Nordbotten].
        double const csp_min = 1.0 / 6.0;
        double const csp_max = 1.0;
        if (std::isnan(fixed_stress_stabilization_parameter) ||
            fixed_stress_stabilization_parameter < csp_min ||
            fixed_stress_stabilization_parameter > csp_max)
        {
            WARN(
                "Value of coupling scheme parameter = {:g} is out of "
                "reasonable range ({:g}, {:g}).",
                fixed_stress_stabilization_parameter, csp_min, csp_max);
        }
    }

    bool const fixed_stress_over_time_step =
        //! \ogs_file_param{prj__processes__process__HYDRO_MECHANICS__coupling_scheme__fixed_stress_over_time_step}
        coupling_scheme_config->getConfigParameter<bool>(
            "fixed_stress_over_time_step", false);

    return Staggered{fixed_stress_stabilization_parameter,
                     fixed_stress_over_time_step};
}

}  // namespace ProcessLib::Common::HydroMechanics
