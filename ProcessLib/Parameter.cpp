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

#include "BaseLib/Error.h"
#include "MeshLib/Elements/Element.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createConstParameter(
    BaseLib::ConfigTree const& config)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "Constant");
    //! \ogs_file_param{parameter__Constant__value}
    auto value = config.getConfigParameter<double>("value");
    DBUG("Using value %g", value);

    return std::unique_ptr<ParameterBase>(new ConstParameter<double>(value));
}

std::unique_ptr<ParameterBase> createMeshPropertyParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{parameter__type}
    config.checkConfigParameter("type", "MeshProperty");
    //! \ogs_file_param{parameter__MeshProperty__field_name}
    auto field_name = config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    if (!mesh.getProperties().hasPropertyVector(field_name))
    {
        OGS_FATAL("The required property %s does not exists in the mesh.",
            field_name.c_str());
    }
    auto const& property =
        mesh.getProperties().template getPropertyVector<double>(field_name);
    if (!property)
    {
        OGS_FATAL("The required property %s is not of the requested type.",
            field_name.c_str());
    }

    return std::unique_ptr<ParameterBase>(
        new MeshElementParameter<double>(*property));
}
}  // namespace ProcessLib
