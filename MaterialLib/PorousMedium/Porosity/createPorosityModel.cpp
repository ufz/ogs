/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createPorosityModel.cpp
 *
 * Created on August 16, 2016, 1:16 PM
 */

#include "createPorosityModel.h"

#include "BaseLib/Error.h"
#include "BaseLib/ConfigTree.h"

#include "ProcessLib/Utils/ProcessUtils.h"

#include "Porosity.h"

namespace MaterialLib
{
namespace PorousMedium
{
std::unique_ptr<Porosity> createPorosityModel(BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{material__porous_medium__porosity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
    {
        auto const& constant_porosity = ProcessLib::findParameter<double>(
            config,
            //! \ogs_file_param_special{material__porous_medium__porosity__porosity_parameter}
            "porosity_parameter", parameters, 1);

        return std::make_unique<Porosity>(constant_porosity);
    }

    OGS_FATAL("The porosity type %s is unavailable.\n",
              "The available type is Constant.",
              type.data());
}

}  // end namespace
}  // end namespace
