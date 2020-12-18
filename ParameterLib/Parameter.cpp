/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "RandomFieldMeshElementParameter.h"
#include "TimeDependentHeterogeneousParameter.h"

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
        INFO("ConstantParameter: {:s}", name);
        return createConstantParameter(name, config);
    }
    if (type == "CurveScaled")
    {
        INFO("CurveScaledParameter: {:s}", name);
        return createCurveScaledParameter(name, config, curves);
    }
    if (type == "Function")
    {
        INFO("FunctionParameter: {:s}", name);
        return createFunctionParameter(name, config, curves);
    }
    if (type == "Group")
    {
        INFO("GroupBasedParameter: {:s}", name);
        return createGroupBasedParameter(name, config, mesh);
    }
    if (type == "MeshElement")
    {
        INFO("MeshElementParameter: {:s}", name);
        return createMeshElementParameter(name, config, mesh);
    }
    if (type == "MeshNode")
    {
        INFO("MeshNodeParameter: {:s}", name);
        return createMeshNodeParameter(name, config, mesh);
    }
    if (type == "RandomFieldMeshElement")
    {
        auto& mesh_var = *BaseLib::findElementOrError(
            begin(meshes), end(meshes),
            [&mesh_name](auto const& m) { return m->getName() == mesh_name; },
            "Expected to find a mesh named " + mesh_name + ".");
        INFO("RandomFieldMeshElement: {:s}", name);
        return createRandomFieldMeshElementParameter(name, config, mesh_var);
    }
    if (type == "TimeDependentHeterogeneousParameter")
    {
        INFO("TimeDependentHeterogeneousParameter: {:s}", name);
        return createTimeDependentHeterogeneousParameter(name, config);
    }

    OGS_FATAL("Cannot construct a parameter of given type '{:s}'.", type);
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
           mesh.getName() +
           "'. The same mesh (the same name) has to be referenced in the "
           "project file. Possible reasons are:\n - the parameter used for the "
           "initial condition is not defined on the bulk mesh,\n - the "
           "parameter's domain of definition mesh differs from the boundary "
           "condition or source term domain of definition mesh.";
}
}  // namespace ParameterLib
