/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshElementParameter.h"

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createMeshElementParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "MeshElement");
    auto const field_name =
        //! \ogs_file_param{prj__parameters__parameter__MeshElement__field_name}
        config.getConfigParameter<std::string>("field_name");
    DBUG("Using field_name {:s}", field_name);

    // TODO other data types than only double
    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Cell)
    {
        OGS_FATAL("The mesh property `{:s}' is not an element property.",
                  field_name);
    }

    return std::make_unique<MeshElementParameter<double>>(name, mesh,
                                                          *property);
}

}  // namespace ParameterLib
