/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   createPermeabilityModel.cpp
 *
 * Created on August 17, 2016, 2:36 PM
 */

#include "createPermeabilityModel.h"

#include <cassert>

#include "BaseLib/ConfigTree.h"

#include "BaseLib/Error.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "ProcessLib/Parameter/ConstantParameter.h"
#include "ProcessLib/Parameter/SpatialPosition.h"
#include "ProcessLib/Utils/ProcessUtils.h"

namespace MaterialLib
{
namespace PorousMedium
{
std::unique_ptr<Permeability> createPermeabilityModel(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{material__porous_medium__permeability__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        auto const& permeability_parameter = ProcessLib::findParameter<double>(
            config,
            //! \ogs_file_param_special{material__porous_medium__permeability__permeability_tensor_entries}
            "permeability_tensor_entries", parameters, 0);

        unsigned dimension = static_cast<unsigned>(
            std::sqrt(permeability_parameter.getNumberOfComponents()));
        if (permeability_parameter.getNumberOfComponents() !=
            dimension * dimension)
        {
            OGS_FATAL(
                "The given parameter has %d components, but the permeability "
                "tensor is defined for a %d dimensional problem.",
                permeability_parameter.getNumberOfComponents(), dimension);
        }

        return std::make_unique<Permeability>(
            permeability_parameter, dimension);
    }
    OGS_FATAL("The permeability type '%s' is unavailable.\n",
              "The available types are \n\tConstant.",
              type.data());
}

}  // end of namespace
}  // end of namespace
