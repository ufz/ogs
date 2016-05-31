/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Parameter.h"

#include <boost/optional.hpp>
#include <logog/include/logog.hpp>

#include "MeshLib/Elements/Element.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createConstParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config)
{
    config.checkConfParam("type", "Constant");
    auto value = config.getConfParam<double>("value");
    DBUG("Using value %g", value);

    return std::unique_ptr<ParameterBase>(
        new ConstParameter<double>(name, value));
}

std::unique_ptr<ParameterBase> createMeshPropertyParameter(
    std::string const& name,
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    config.checkConfParam("type", "MeshProperty");
    auto field_name = config.getConfParam<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    if (!mesh.getProperties().hasPropertyVector(field_name))
    {
        ERR("The required property %s does not exists in the mesh.",
            field_name.c_str());
        std::abort();
    }
    auto const& property =
        mesh.getProperties().template getPropertyVector<double>(field_name);
    if (!property)
    {
        ERR("The required property %s is not of the requested type.",
            field_name.c_str());
        std::abort();
    }

    return std::unique_ptr<ParameterBase>(
        new MeshPropertyParameter<double>(name, *property));
}
}  // namespace ProcessLib
