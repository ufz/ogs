/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file:   createStorageModel.cpp
 *
 * Created on August 16, 2016, 1:16 PM
 */

#include "createStorageModel.h"

#include "BaseLib/Error.h"

#include "ConstantStorage.h"

namespace MaterialLib
{
namespace PorousMedium
{
Storage* createStorageModel(BaseLib::ConfigTree const* const config)
{
    //! \ogs_file_param{material__porous_medium__storage__type}
    auto const type = config->getConfigParameter<std::string>("type");

    if (type.find("constant") != std::string::npos)
        //! \ogs_file_param{material__fluid__viscosity__value}
        return new ConstantStorage(config->getConfigParameter<double>("value"));
    else
    {
        OGS_FATAL(
            "The storage type is unavailable.\n"
            "The available type is constant.");
    }
}

}  // end namespace
}  // end namespace
