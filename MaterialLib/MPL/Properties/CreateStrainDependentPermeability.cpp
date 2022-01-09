/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 10:08 AM
 */

#include "CreateStrainDependentPermeability.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "Parameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "ParameterLib/Parameter.h"
#include "ParameterLib/Utils.h"
#include "StrainDependentPermeability.h"

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createStrainDependentPermeability(
    int const geometry_dimension,
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system)
{
    if ((geometry_dimension != 2) && (geometry_dimension != 3))
    {
        OGS_FATAL(
            "The StrainDependentPermeability is implemented only for 2D or 3D "
            "problems");
    }

    //! \ogs_file_param{properties__property__type}
    config.checkConfigParameter("type", "StrainDependentPermeability");

    // Second access for storage.
    //! \ogs_file_param{properties__property__name}
    auto property_name = config.peekConfigParameter<std::string>("name");

    DBUG("Create StrainDependentPermeability property {:s}.", property_name);

    std::string const& parameter_name =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__initial_permeability}
        config.getConfigParameter<std::string>("initial_permeability");
    auto const& parameter_k0 = ParameterLib::findParameter<double>(
        parameter_name, parameters, 0, nullptr);

    auto const b1 =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__b1}
        config.getConfigParameter<double>("b1");
    auto const b2 =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__b2}
        config.getConfigParameter<double>("b2");
    auto const b3 =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__b3}
        config.getConfigParameter<double>("b3");
    auto const minimum_permeability =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__minimum_permeability}
        config.getConfigParameter<double>("minimum_permeability");
    auto const maximum_permeability =
        //! \ogs_file_param{properties__property__StrainDependentPermeability__maximum_permeability}
        config.getConfigParameter<double>("maximum_permeability");

    if (minimum_permeability > maximum_permeability)
    {
        OGS_FATAL(
            "The value of minimum_permeability of {:e} is larger that the "
            "value of maximum_permeability of {:e} in "
            "StrainDependentPermeability",
            minimum_permeability, maximum_permeability);
    }

    if (geometry_dimension == 2)
    {
        return std::make_unique<StrainDependentPermeability<2>>(
            std::move(property_name), parameter_k0, b1, b2, b3,
            minimum_permeability, maximum_permeability,
            local_coordinate_system);
    }

    return std::make_unique<StrainDependentPermeability<3>>(
        std::move(property_name), parameter_k0, b1, b2, b3,
        minimum_permeability, maximum_permeability, local_coordinate_system);
}
}  // namespace MaterialPropertyLib
