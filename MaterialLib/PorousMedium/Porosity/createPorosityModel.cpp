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
#include "BaseLib/ConfigTree.h"

#include "Porosity.h"
#include "ConstantPorosity.h"

namespace MaterialLib
{
namespace PorousMedium
{
std::unique_ptr<Porosity> createPorosityModel(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__porosity__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::unique_ptr<Porosity>(
            //! \ogs_file_param{material__porous_medium__porosity__Constant__value}
            new ConstantPorosity(config.getConfigParameter<double>("value")));
    else
    {
        OGS_FATAL("The porosity type %s is unavailable.\n",
                  "The available type is \n\tConstant.", type.data());
    }
}

}  // end namespace
}  // end namespace
