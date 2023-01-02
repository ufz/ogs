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
NumericalStabilization::~NumericalStabilization() {}

IsotropicDiffusionStabilization::IsotropicDiffusionStabilization(
    double const cutoff_velocity,
    double const tuning_parameter,
    std::vector<double>&& element_sizes_vector)
    : NumericalStabilization(cutoff_velocity),
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
    std::size_t const elemend_id, double const velocity_norm) const
{
    if (velocity_norm < cutoff_velocity_)
    {
        return 0.0;
    }
    return 0.5 * tuning_parameter_ * velocity_norm * element_sizes_[elemend_id];
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

    const auto type =
        //! \ogs_file_param{prj__processes__process__numerical_stabilization__type}
        stabilization_config->getConfigParameter<std::string>("type");
    if (type.compare("IsotropicDiffusion") == 0)
    {
        INFO("Using numerical stabilization {:s}.", type);

        auto const stabilization_cutoff_velocity_opt =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__IsotropicDiffusion__cutoff_velocity}
            stabilization_config->getConfigParameterOptional<double>(
                "cutoff_velocity");

        const auto tuning_parameter =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__IsotropicDiffusion__tuning_parameter}
            stabilization_config->getConfigParameter<double>(
                "tuning_parameter");

        return std::make_unique<IsotropicDiffusionStabilization>(
            *stabilization_cutoff_velocity_opt,
            tuning_parameter,
            MeshLib::getMaxiumElementEdgeLengths(mesh.getElements()));
    }
    else if (type.compare("FullUpwind") == 0)
    {
        INFO("Numerical stabilization of {:s} is used", type);

        auto const stabilization_cutoff_velocity_ptr =
            //! \ogs_file_param{prj__processes__process__numerical_stabilization__FullUpwind__cutoff_velocity}
            stabilization_config->getConfigParameterOptional<double>(
                "cutoff_velocity");

        return std::make_unique<FullUpwind>(*stabilization_cutoff_velocity_ptr);
    }

    OGS_FATAL("The stabilization type {:s} is not available.", type);
}

}  // namespace NumLib
