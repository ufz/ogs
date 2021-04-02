/**
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on August 17, 2016, 2:36 PM
 */

#include "createPermeabilityModel.h"

#include <cassert>

#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MaterialLib/PorousMedium/Permeability/DupuitPermeability.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/SpatialPosition.h"
#include "ParameterLib/Utils.h"

namespace MaterialLib
{
namespace PorousMedium
{
std::unique_ptr<Permeability> createPermeabilityModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{material__porous_medium__permeability__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        auto const& permeability_parameter = ParameterLib::findParameter<
            double>(
            config,
            //! \ogs_file_param_special{material__porous_medium__permeability__permeability_tensor_entries}
            "permeability_tensor_entries", parameters, 0);

        int dimension = static_cast<int>(
            std::sqrt(permeability_parameter.getNumberOfGlobalComponents()));
        if (permeability_parameter.getNumberOfGlobalComponents() !=
            dimension * dimension)
        {
            OGS_FATAL(
                "The given parameter has {:d} components, but the permeability "
                "tensor is defined for a {:d} dimensional problem.",
                permeability_parameter.getNumberOfGlobalComponents(),
                dimension);
        }

        return std::make_unique<Permeability>(permeability_parameter,
                                              dimension);
    }

    if (type == "Dupuit")
    {
        auto const& permeability_parameter = ParameterLib::findParameter<
            double>(
            config,
            //! \ogs_file_param_special{material__porous_medium__permeability__permeability_tensor_entries}
            "permeability_tensor_entries", parameters, 0);

        int dimension = static_cast<int>(
            std::sqrt(permeability_parameter.getNumberOfGlobalComponents()));
        if (permeability_parameter.getNumberOfGlobalComponents() !=
            dimension * dimension)
        {
            OGS_FATAL(
                "The given parameter has {:d} components, but the permeability "
                "tensor is defined for a {:d} dimensional problem.",
                permeability_parameter.getNumberOfGlobalComponents(),
                dimension);
        }

        return std::make_unique<DupuitPermeability>(permeability_parameter,
                                                    dimension);
    }
    OGS_FATAL("The permeability type '{:s}' is unavailable.\n",
              "The available types are \n\tConstant.",
              type.data());
}

}  // namespace PorousMedium
}  // namespace MaterialLib
