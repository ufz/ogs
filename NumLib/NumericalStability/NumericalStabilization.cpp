/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 2, 2022, 2:40 PM
 */
#include "NumericalStabilization.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/getMaxiumElementEdgeLengths.h"

namespace NumLib
{

IsotropicDiffusionStabilization::IsotropicDiffusionStabilization(
    double const cutoff_velocity,
    double const tuning_parameter,
    std::vector<double>&& element_sizes_vector)
    : cutoff_velocity_(cutoff_velocity),
      tuning_parameter_(tuning_parameter),
      element_sizes_(std::move(element_sizes_vector))
{
    if (tuning_parameter_ < 0 || tuning_parameter_ > 1.0)
    {
        OGS_FATAL(
            "The tuning parameter value {:g} for "
            "IsotropicDiffusion stabilization is out of range [0, 1]",
            tuning_parameter_);
    }
}

double IsotropicDiffusionStabilization::computeArtificialDiffusion(
    std::size_t const element_id, double const velocity_norm) const
{
    if (velocity_norm < cutoff_velocity_)
    {
        return 0.0;
    }
    return 0.5 * tuning_parameter_ * velocity_norm * element_sizes_[element_id];
}

std::unique_ptr<NumericalStabilization> createNumericalStabilization(
    MeshLib::Mesh const& mesh, BaseLib::ConfigTree const& config)
{
    auto const stabilization_config =
        //! \ogs_file_param{prj__processes__process__numerical_stabilization}
        config.getConfigSubtreeOptional("numerical_stabilization");
    if (!stabilization_config)
    {
        return nullptr;
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

        return std::make_unique<IsotropicDiffusionStabilization>(
            cutoff_velocity,
            tuning_parameter,
            MeshLib::getMaxiumElementEdgeLengths(mesh.getElements()));
    }
    if (type == "FullUpwind")
    {
        auto const cutoff_velocity =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__FullUpwind__cutoff_velocity}
            stabilization_config->getConfigParameter<double>("cutoff_velocity");

        return std::make_unique<FullUpwind>(cutoff_velocity);
    }

    OGS_FATAL("The stabilization type {:s} is not available.", type);
}

}  // namespace NumLib
