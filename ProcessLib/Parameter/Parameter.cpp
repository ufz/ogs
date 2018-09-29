/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Parameter.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"

#include "ConstantParameter.h"
#include "CurveScaledParameter.h"
#include "GroupBasedParameter.h"
#include "MeshElementParameter.h"
#include "MeshNodeParameter.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createParameter(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& meshes,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves)
{

    //! \ogs_file_param{prj__parameters__parameter__name}
    auto const name = config.getConfigParameter<std::string>("name");
    //! \ogs_file_param{prj__parameters__parameter__type}
    auto const type = config.peekConfigParameter<std::string>("type");

    // Create parameter based on the provided type.
    if (type == "Constant")
    {
        INFO("ConstantParameter: %s", name.c_str());
        auto param = createConstantParameter(name, config);
        return param;
    }
    if (type == "CurveScaled")
    {
        INFO("CurveScaledParameter: %s", name.c_str());
        auto param = createCurveScaledParameter(name, config, curves);
        return param;
    }
    if (type == "Group")
    {
        INFO("GroupBasedParameter: %s", name.c_str());
        auto param = createGroupBasedParameter(name, config, *meshes.front());
        return param;
    }
    if (type == "MeshElement")
    {
        INFO("MeshElementParameter: %s", name.c_str());
        auto param = createMeshElementParameter(name, config, *meshes.front());
        return param;
    }
    if (type == "MeshNode")
    {
        INFO("MeshNodeParameter: %s", name.c_str());
        auto param = createMeshNodeParameter(name, config, *meshes.front());
        return param;
    }

    OGS_FATAL("Cannot construct a parameter of given type \'%s\'.",
              type.c_str());
}

}  // ProcessLib
