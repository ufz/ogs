/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SourceTermBuilder.h"
#include "SourceTermConfig.h"
#include "CreateNodalSourceTerm.h"
#include "NodalSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<NodalSourceTerm> SourceTermBuilder::createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned integration_order,
    const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    if (type == "Nodal")
    {
        return createNodalSourceTerm(config, dof_table, mesh, variable_id,
                                     integration_order, shapefunction_order,
                                     parameters);
    }

    OGS_FATAL("Unknown source term type: `%s'.", type.c_str());
}

std::unique_ptr<NodalSourceTerm> SourceTermBuilder::createNodalSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, const MeshLib::Mesh& mesh,
    const int variable_id, const unsigned /*integration_order*/,
    const unsigned /*shapefunction_order*/,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    return ProcessLib::createNodalSourceTerm(
        config.config, config.mesh, dof_table, mesh.getID(), variable_id,
        *config.component_id, parameters);
}

}  // ProcessLib
