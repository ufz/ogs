/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateNodalSourceTerm.h"

#include "BaseLib/FileTools.h"
#include "ProcessLib/Utils/ProcessUtils.h"
#include "NodalSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<NodalSourceTerm> createNodalSourceTerm(
    BaseLib::ConfigTree const& config,
    const NumLib::LocalToGlobalIndexMap& dof_table, std::size_t const mesh_id,
    std::size_t const node_id, const int variable_id, const int component_id)
{
    DBUG("Constructing NodalSourceTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Nodal");

    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Nodal__value}
    auto const nodal_value = config.getConfigParameter<double>("value");
    DBUG("Using value %f as nodal source term", nodal_value);

    return std::make_unique<NodalSourceTerm>(
        dof_table, mesh_id, node_id, variable_id, component_id, nodal_value);
}

}  // namespace ProcessLib
