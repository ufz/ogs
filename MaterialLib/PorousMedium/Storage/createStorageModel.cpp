/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "BaseLib/ConfigTree.h"

#include "Storage.h"
#include "ConstantStorage.h"

namespace MaterialLib
{
namespace PorousMedium
{
std::unique_ptr<Storage> createStorageModel(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__storage__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "Constant")
        return std::make_unique<ConstantStorage>(
            //! \ogs_file_param{material__porous_medium__storage__Constant__value}
            config.getConfigParameter<double>("value"));

    OGS_FATAL("The storage type %s is unavailable.\n", type.data(),
              "The available type is Constant.");
}

}  // end namespace
}  // end namespace
