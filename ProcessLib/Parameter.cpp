/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Parameter.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ConstantParameter.h"
#include "MeshElementParameter.h"
#include "MeshNodeParameter.h"

namespace ProcessLib
{

std::unique_ptr<ParameterBase> createParameter(
    BaseLib::ConfigTree const& config, std::vector<MeshLib::Mesh*> const& meshes)
{

    //! \ogs_file_param{parameter__name}
    auto name = config.getConfigParameter<std::string>("name");
    //! \ogs_file_param{parameter__type}
    auto type = config.peekConfigParameter<std::string>("type");

    // Create parameter based on the provided type.
    if (type == "Constant")
    {
        INFO("ConstantParameter: %s", name.c_str());
        auto param = createConstantParameter(config);
        param->name = name;
        return param;
    }
    else if (type == "MeshElement")
    {
        INFO("MeshElementParameter: %s", name.c_str());
        auto param = createMeshElementParameter(config, *meshes.front());
        param->name = name;
        return param;
    }
    else if (type == "MeshNode")
    {
        INFO("MeshElementParameter: %s", name.c_str());
        auto param = createMeshNodeParameter(config, *meshes.front());
        param->name = name;
        return param;
    }
    else
    {
        OGS_FATAL("Cannot construct property of given type \'%s\'.",
                  type.c_str());
    }
}

}  // ProcessLib
