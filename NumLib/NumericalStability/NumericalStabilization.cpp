/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::vector<double>&& element_sizes)
    : cutoff_velocity_(cutoff_velocity),
      tuning_parameter_(tuning_parameter),
      element_sizes_(std::move(element_sizes))
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
}  // namespace NumLib
