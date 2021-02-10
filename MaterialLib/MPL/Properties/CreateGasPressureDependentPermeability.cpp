/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 12, 2021, x:xx AM
 */

#include "CreateGasPressureDependentPermeability.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "GasPressureDependentPermeability.h"
#include "Parameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createGasPressureDependentPermeability(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    if ((geometry_dimension != 2) && (geometry_dimension != 3))
    {
        OGS_FATAL(
            "The GasPressureDependentPermeability is implemented only "
            "for 2D or 3D problems");
    }

    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "GasPressureDependentPermeability");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create GasPressureDependentPermeability property {:s}.",
         property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__initial_permeability}
        config.getConfigParameter<std::string>("initial_permeability");
    auto const& parameter_k0 = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);

    auto const a1 =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__b1}
        config.getConfigParameter<double>("a1");
    auto const a2 =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__b2}
        config.getConfigParameter<double>("a2");
    auto const pthr =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__b3}
        config.getConfigParameter<double>("pthr");
    auto const minimum_permeability =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__minimum_permeability}
        config.getConfigParameter<double>("minimum_permeability");
    auto const maximum_permeability =
        //! \ogs_file_param{properties__property__GasPressureDependentPermeability__maximum_permeability}
        config.getConfigParameter<double>("maximum_permeability");

    if (minimum_permeability > maximum_permeability)
    {
        OGS_FATAL(
            "The value of minimum_permeability of {:e} is larger that the "
            "value of maximum_permeability of {:e} in "
            "GasPressureDependentPermeability",
            minimum_permeability, maximum_permeability);
    }

    if (geometry_dimension == 2)
    {
        return std::make_unique<GasPressureDependentPermeability<2>>(
            std::move(property_name), parameter_k0, a1, a2, pthr,
            minimum_permeability, maximum_permeability,
            local_coordinate_system);
    }

    return std::make_unique<GasPressureDependentPermeability<3>>(
        std::move(property_name), parameter_k0, a1, a2, pthr,
        minimum_permeability, maximum_permeability, local_coordinate_system);
}
}  // namespace MaterialPropertyLib
