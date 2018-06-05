/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
    BaseLib::ConfigTree const& config, MeshLib::Mesh& st_mesh,
    const NumLib::LocalToGlobalIndexMap& dof_table,
    std::size_t const bulk_mesh_id, const int variable_id,
    const int component_id,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing NodalSourceTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Nodal");

    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Nodal__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter %s as nodal source term.", param_name.c_str());

    auto& param = findParameter<double>(param_name, parameters, 1);

    return std::make_unique<NodalSourceTerm>(dof_table, bulk_mesh_id, st_mesh,
                                             variable_id, component_id, param);
}

}  // namespace ProcessLib
