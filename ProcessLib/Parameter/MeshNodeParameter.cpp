/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshNodeParameter.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"

namespace ProcessLib
{
std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "MeshNode");
    //! \ogs_file_param{prj__parameters__parameter__MeshNode__field_name}
    auto const field_name = config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    // TODO other data types than only double
    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node) {
        OGS_FATAL("The mesh property `%s' is not a nodal property.",
                  field_name.c_str());
    }

    return std::unique_ptr<ParameterBase>(
        new MeshNodeParameter<double>(name, *property));
}

}  // ProcessLib
