/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "ConstantPorosity.h"

namespace MaterialLib
{
namespace PorousMedium
{
Porosity* createPorosityModel(BaseLib::ConfigTree const* const config)
{
    //! \ogs_file_param{material__porous_medium__storage__type}
    auto const type = config->getConfigParameter<std::string>("type");

    if (type.find("constant") != std::string::npos)
        //! \ogs_file_param{material__porous_medium__porosity__value}
        return new ConstantPorosity(
            config->getConfigParameter<double>("value"));
    else
    {
        OGS_FATAL(
            "The storage type is unavailable.\n"
            "The available type is constant.");
    }
}

}  // end namespace
}  // end namespace
