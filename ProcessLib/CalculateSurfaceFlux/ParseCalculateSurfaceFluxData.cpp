/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ParseCalculateSurfaceFluxData.h"

namespace ProcessLib
{
void parseCalculateSurfaceFluxData(BaseLib::ConfigTree const& config,
                                   std::string& mesh_name,
                                   std::string& property_name,
                                   std::string& out_fname)
{
    auto calculatesurfaceflux_config =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux}
        config.getConfigSubtreeOptional("calculatesurfaceflux");
    if (!calculatesurfaceflux_config)
        return;

    mesh_name =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__mesh}
        calculatesurfaceflux_config->getConfigParameter<std::string>("mesh");
    property_name =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__property_name}
        calculatesurfaceflux_config->getConfigParameter<std::string>(
            "property_name");
    out_fname =
        //! \ogs_file_param{prj__processes__process__calculatesurfaceflux__output_mesh}
        calculatesurfaceflux_config->getConfigParameter<std::string>(
            "output_mesh");
}

}  // namespace ProcessLib
