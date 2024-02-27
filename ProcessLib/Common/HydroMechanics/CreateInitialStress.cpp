/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 13, 2024, 10:20 AM
 */

#include "CreateInitialStress.h"

#include <string>

#include "BaseLib/ConfigTree.h"
#include "InitialStress.h"
#include "MathLib/KelvinVector.h"
#include "MeshLib/Mesh.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
template <int DisplacementDim>
InitialStress createInitialStress(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type)
{
    auto const config_stress0 =
        //! \ogs_file_param{prj__processes__process__initial_stress}
        config.getConfigSubtreeOptional("initial_stress");

    if (!config_stress0)
    {
        return {};
    }

    auto const stress0_type_opt =
        //! \ogs_file_attr{prj__processes__process__initial_stress__type}
        config_stress0->getConfigAttributeOptional<std::string>("type");

    // If the type is mandatory required but it is not presented:
    if (mandatory_stress_type && !stress0_type_opt)
    {
        OGS_FATAL(
            "The attribute \"type\" of the tag \"initial_stress\" is required "
            "for the current process setting, for example using subtype = "
            "\"StressSaturation_StrainPressureTemperature\" in TRM. "
            "However it is missing.");
    }

    auto const stress0_type_str =
        stress0_type_opt ? *stress0_type_opt : "effective";

    InitialStress::Type stress0_type;
    if (stress0_type_str == "total")
    {
        stress0_type = InitialStress::Type::Total;
    }
    else if (stress0_type_str == "effective")
    {
        stress0_type = InitialStress::Type::Effective;
    }
    else
    {
        OGS_FATAL(
            "The initial stress type must be \"total\" or "
            "\"effective\". But the given one is {:s}",
            stress0_type_str);
    }

    auto const parameter_name = config_stress0->getValue<std::string>();
    auto const initial_stress = &ParameterLib::findParameter<double>(
        parameter_name, parameters,
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim),
        &mesh);

    return {initial_stress, stress0_type};
}

template InitialStress createInitialStress<2>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type);
template InitialStress createInitialStress<3>(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::Mesh const& mesh, bool const mandatory_stress_type);

};  // namespace ProcessLib
