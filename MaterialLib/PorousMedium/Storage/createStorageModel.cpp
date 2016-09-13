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
std::unique_ptr<Storage> createStorageModel(BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{material__porous_medium__storage__type}
    auto const type = config.getConfigParameter<std::string>("type");

    if (type == "constant")
        //! \ogs_file_param{material__porous_medium__storage__constant__value}
        return std::unique_ptr<Storage>(
                 new ConstantStorage(config.getConfigParameter<double>("value")) );
    else
    {
        OGS_FATAL(
            "The storage type %s is unavailable.\n", type.data(),
            "The available type is constant.");
    }
}

}  // end namespace
}  // end namespace
