/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "ProcessLib/Utils/ProcessUtils.h"
#include "VolumetricSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createVolumetricSourceTerm(
    BaseLib::ConfigTree const& config,
    MeshLib::Mesh const& source_term_mesh,
    NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order, unsigned const shapefunction_order,
    int const variable_id, int const component_id)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Volumetric");

    DBUG("Constructing VolumetricSourceTerm from config.");

    // source term field name
    auto const& volumetric_source_term_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Volumetric__parameter}
        config.getConfigParameter<std::string>("parameter");
    auto& volumetric_source_term = findParameter<double>(
        volumetric_source_term_parameter_name, parameters, 1);

    DBUG("Using '%s` as volumetric source term parameter.",
         volumetric_source_term.name.c_str());

    return std::make_unique<VolumetricSourceTerm>(
        source_term_mesh, source_term_dof_table, integration_order,
        shapefunction_order, variable_id, component_id,
        volumetric_source_term);
}

}  // namespace ProcessLib
