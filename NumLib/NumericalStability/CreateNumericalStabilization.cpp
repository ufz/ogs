/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CreateNumericalStabilization.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getMaxiumElementEdgeLengths.h"
#include "NumericalStabilization.h"

namespace NumLib
{
NumericalStabilization createNumericalStabilization(
    MeshLib::Mesh const& mesh, BaseLib::ConfigTree const& config)
{
    auto const stabilization_config =
        //! \ogs_file_param{prj__processes__process__numerical_stabilization}
        config.getConfigSubtreeOptional("numerical_stabilization");
    if (!stabilization_config)
    {
        return NoStabilization{};
    }

    auto const type =
        //! \ogs_file_param{prj__processes__process__numerical_stabilization__type}
        stabilization_config->getConfigParameter<std::string>("type");

    INFO("Using {:s} numerical stabilization.", type);
    if (type == "IsotropicDiffusion")
    {
        auto const cutoff_velocity =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__IsotropicDiffusion__cutoff_velocity}
            stabilization_config->getConfigParameter<double>("cutoff_velocity");

        auto const tuning_parameter =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__IsotropicDiffusion__tuning_parameter}
            stabilization_config->getConfigParameter<double>(
                "tuning_parameter");

        return IsotropicDiffusionStabilization{
            cutoff_velocity, tuning_parameter,
            MeshLib::getMaxiumElementEdgeLengths(mesh.getElements())};
    }
    if (type == "FullUpwind")
    {
        auto const cutoff_velocity =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__FullUpwind__cutoff_velocity}
            stabilization_config->getConfigParameter<double>("cutoff_velocity");

        return FullUpwind{cutoff_velocity};
    }
    if (type == "FluxCorrectedTransport")
    {
        return FluxCorrectedTransport();
    }

    OGS_FATAL("The stabilization type {:s} is not available.", type);
}
}  // namespace NumLib
