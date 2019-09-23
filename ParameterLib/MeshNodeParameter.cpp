/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshNodeParameter.h"
#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createMeshNodeParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "MeshNode");
    auto const field_name =
        //! \ogs_file_param{prj__parameters__parameter__MeshNode__field_name}
        config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name %s", field_name.c_str());

    // TODO other data types than only double
    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Node)
    {
        OGS_FATAL("The mesh property `%s' is not a nodal property.",
                  field_name.c_str());
    }

    return std::make_unique<MeshNodeParameter<double>>(name, mesh, *property);
}

}  // namespace ParameterLib
