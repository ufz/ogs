/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateVolumetricSourceTerm.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Utils.h"
#include "VolumetricSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createVolumetricSourceTerm(
    BaseLib::ConfigTree const& config, unsigned const bulk_mesh_dimension,
    MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order, unsigned const shapefunction_order)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.peekConfigParameter<std::string>("type");
    if (type == "Line")
    {
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
        config.checkConfigParameter("type", "Line");
        DBUG("Constructing LineSourceTerm from config.");
    }
    else
    {
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
        config.checkConfigParameter("type", "Volumetric");
        DBUG("Constructing VolumetricSourceTerm from config.");
    }

    // source term field name
    auto const& volumetric_source_term_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Volumetric__parameter}
        config.getConfigParameter<std::string>("parameter");
    auto& volumetric_source_term = ParameterLib::findParameter<double>(
        volumetric_source_term_parameter_name, parameters, 1,
        &source_term_mesh);

    DBUG("Using '{:s}' as volumetric source term parameter.",
         volumetric_source_term.name);

    return std::make_unique<VolumetricSourceTerm>(
        bulk_mesh_dimension, source_term_mesh, std::move(source_term_dof_table),
        integration_order, shapefunction_order, volumetric_source_term);
}

}  // namespace ProcessLib
