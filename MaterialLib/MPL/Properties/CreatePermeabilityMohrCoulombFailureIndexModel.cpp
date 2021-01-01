/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 5, 2020, 9:07 AM
 */

#include "CreatePermeabilityMohrCoulombFailureIndexModel.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "Parameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "PermeabilityMohrCoulombFailureIndexModel.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createPermeabilityMohrCoulombFailureIndexModel(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    if ((geometry_dimension != 2) && (geometry_dimension != 3))
    {
        OGS_FATAL(
            "The PermeabilityMohrCoulombFailureIndexModel is implemented only "
            "for 2D or 3D problems");
    }

    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type",
                                "PermeabilityMohrCoulombFailureIndexModel");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create PermeabilityMohrCoulombFailureIndexModel property {:s}.",
         property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__initial_permeability}
        config.getConfigParameter<std::string>(
            "initial_permeability");
    auto const& parameter_k0 = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);

    auto const kr =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__reference_permeability}
        config.getConfigParameter<double>("reference_permeability");
    auto const b =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__fitting_factor}
        config.getConfigParameter<double>("fitting_factor");
    auto const c =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__cohesion}
        config.getConfigParameter<double>("cohesion");
    auto const phi =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__friction_angle}
        config.getConfigParameter<double>("friction_angle");
    auto const max_k =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__maximum_permeability}
        config.getConfigParameter<double>("maximum_permeability");
    auto const t_sigma_max =
        //! \ogs_file_param{properties__property__PermeabilityMohrCoulombFailureIndexModel__tensile_strength_parameter}
        config.getConfigParameter<double>("tensile_strength_parameter");

    if (geometry_dimension == 2)
    {
        return std::make_unique<PermeabilityMohrCoulombFailureIndexModel<2>>(
            std::move(property_name), parameter_k0, kr, b, c, phi, max_k,
            t_sigma_max, local_coordinate_system);
    }

    return std::make_unique<PermeabilityMohrCoulombFailureIndexModel<3>>(
        std::move(property_name), parameter_k0, kr, b, c, phi, max_k,
        t_sigma_max, local_coordinate_system);
}
}  // namespace MaterialPropertyLib
