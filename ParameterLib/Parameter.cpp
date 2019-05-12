/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "FunctionParameter.h"
#include "GroupBasedParameter.h"
#include "MeshElementParameter.h"
#include "MeshNodeParameter.h"

namespace ParameterLib
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

    // Either the mesh name is given, or the first mesh's name will be
    // taken.
    auto const mesh_name =
        //! \ogs_file_param{prj__parameters__parameter__mesh}
        config.getConfigParameter<std::string>("mesh", meshes[0]->getName());

    auto const& mesh = *BaseLib::findElementOrError(
        begin(meshes), end(meshes),
        [&mesh_name](auto const& m) { return m->getName() == mesh_name; },
        "Expected to find a mesh named " + mesh_name + ".");

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
    if (type == "Function")
    {
        INFO("FunctionParameter: %s", name.c_str());
        auto param = createFunctionParameter(name, config, mesh);
        return param;
    }
    if (type == "Group")
    {
        INFO("GroupBasedParameter: %s", name.c_str());
        auto param = createGroupBasedParameter(name, config, mesh);
        return param;
    }
    if (type == "MeshElement")
    {
        INFO("MeshElementParameter: %s", name.c_str());
        auto param = createMeshElementParameter(name, config, mesh);
        return param;
    }
    if (type == "MeshNode")
    {
        INFO("MeshNodeParameter: %s", name.c_str());
        auto param = createMeshNodeParameter(name, config, mesh);
        return param;
    }

    OGS_FATAL("Cannot construct a parameter of given type '%s'.", type.c_str());
}

boost::optional<std::string> isDefinedOnSameMesh(ParameterBase const& parameter,
                                                 MeshLib::Mesh const& mesh)
{
    // Arbitrary domain of definition.
    if (parameter.mesh() == nullptr)
    {
        return {};
    }

    // Equal meshes.
    if (*parameter.mesh() == mesh)
    {
        return {};
    }

    return "The parameter's domain of definition mesh '" +
           parameter.mesh()->getName() + "' differs from the used mesh '" +
           mesh.getName() + "'. Both meshes must be equal.";
}
}  // namespace ParameterLib
