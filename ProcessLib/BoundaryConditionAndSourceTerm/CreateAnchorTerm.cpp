/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateAnchorTerm.h"

#include "AnchorTerm.h"
#include "BaseLib/ConfigTree.h"
#include "BaseLib/Logging.h"
#include "ParameterLib/Utils.h"

namespace ProcessLib
{
template <int GlobalDim>
std::unique_ptr<SourceTerm> createAnchorTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters)
{
    DBUG("Constructing AnchorTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Anchor");

    auto const param_name =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Anchor__anchor_force_constant}
        config.getConfigParameter<std::string>("anchor_force_constant");
    DBUG("Using parameter {:s} as anchor stiffness constant.", param_name);

    auto& param = ParameterLib::findParameter<double>(param_name, parameters, 1,
                                                      &st_mesh);

    for (MeshLib::Element const* const element : st_mesh.getElements())
    {
        if (element->getNumberOfNodes() != 2)
        {
            OGS_FATAL(
                "Every anchor element needs to have precisely two nodes.");
        }
    }

    return std::make_unique<AnchorTerm<GlobalDim>>(
        std::move(dof_table), source_term_mesh_id, st_mesh, variable_id, param);
}

template std::unique_ptr<SourceTerm> createAnchorTerm<2>(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

template std::unique_ptr<SourceTerm> createAnchorTerm<3>(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& st_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    std::size_t const source_term_mesh_id, const int variable_id,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace ProcessLib
