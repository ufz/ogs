// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateSigmoid.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Properties/Sigmoid.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Sigmoid> createSigmoid(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "Sigmoid");

    //! \ogs_file_param{properties__property__name}
    auto const name = config.peekConfigParameter<std::string>("name");
    DBUG("Create sigmoid property '{:s}'.", name);

    auto const steepness =
        //! \ogs_file_param{properties__property__Sigmoid__steepness}
        config.getConfigParameter<double>("steepness");

    auto const midpoint =
        //! \ogs_file_param{properties__property__Sigmoid__midpoint}
        config.getConfigParameter<double>("midpoint");

    auto const lower_bound =
        //! \ogs_file_param{properties__property__Sigmoid__lower_bound}
        config.getConfigParameter<double>("lower_bound");

    auto const upper_bound =
        //! \ogs_file_param{properties__property__Sigmoid__upper_bound}
        config.getConfigParameter<double>("upper_bound");

    // Parse independent variable
    auto const independent_variable_str =
        //! \ogs_file_param{properties__property__Sigmoid__independent_variable}
        config.getConfigParameter<std::string>("independent_variable");

    // Spatial/temporal coordinates not yet supported; would need to
    // reintroduce a SpatialTemporalVariable and dispatch in value().
    if (std::ranges::any_of(std::array{"t", "x", "y", "z"},
                            [&](const char* v)
                            { return independent_variable_str == v; }))

    {
        OGS_FATAL(
            "In the sigmoid property '{:s}', the spatial/temporal independent "
            "variable '{:s}' is not supported yet. Use a material variable "
            "instead. See CreateSigmoid.cpp for how to re-enable "
            "spatial/temporal independent variables.",
            name, independent_variable_str);
    }

    auto const independent_variable =
        convertStringToVariable(independent_variable_str);

    return std::make_unique<Sigmoid>(name, steepness, midpoint, lower_bound,
                                     upper_bound, independent_variable);
}

}  // namespace MaterialPropertyLib
