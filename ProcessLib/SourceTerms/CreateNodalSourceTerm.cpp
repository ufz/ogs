/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateNodalSourceTerm.h"

#include "BaseLib/Logging.h"

#include "BaseLib/ConfigTree.h"
#include "ParameterLib/Utils.h"

#include "NodalSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createNodalSourceTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    const int component_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing NodalSourceTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Nodal");

    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Nodal__parameter}
    auto const param_name = config.getConfigParameter<std::string>("parameter");
    DBUG("Using parameter {:s} as nodal source term.", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &st_mesh);

    return std::make_unique<NodalSourceTerm>(std::move(dof_table),
                                             source_term_mesh_id, st_mesh,
                                             variable_id, component_id, param);
}

}  // namespace ProcessLib
