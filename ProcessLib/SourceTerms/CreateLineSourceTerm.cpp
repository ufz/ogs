/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateLineSourceTerm.h"

#include "BaseLib/ConfigTree.h"
#include "BaseLib/FileTools.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ParameterLib/Utils.h"
#include "LineSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createLineSourceTerm(
    BaseLib::ConfigTree const& config, unsigned const bulk_mesh_dimension,
    MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order, unsigned const shapefunction_order)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Line");

    DBUG("Constructing LineSourceTerm from config.");

    // source term field name
    auto const& line_source_term_parameter_name =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Line__parameter}
        config.getConfigParameter<std::string>("parameter");
    auto& line_source_term = ParameterLib::findParameter<double>(
        line_source_term_parameter_name, parameters, 1, &source_term_mesh);

    DBUG("Using '{:s}' as line source term parameter.",
         line_source_term.name.c_str());

    return std::make_unique<LineSourceTerm>(
        bulk_mesh_dimension, source_term_mesh, std::move(source_term_dof_table),
        integration_order, shapefunction_order, line_source_term);
}

}  // namespace ProcessLib
