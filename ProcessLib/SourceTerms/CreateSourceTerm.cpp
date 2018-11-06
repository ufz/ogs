/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSourceTerm.h"

#include "CreateNodalSourceTerm.h"
#include "CreateVolumetricSourceTerm.h"
#ifdef OGS_USE_PYTHON
#include "Python/CreatePythonSourceTerm.h"
#endif
#include "SourceTerm.h"
#include "SourceTermConfig.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table_bulk,
    const MeshLib::Mesh& source_term_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");


    MeshLib::MeshSubset source_term_mesh_subset(source_term_mesh,
                                                source_term_nodes);


    if (type == "Nodal")
    {
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_source_term(
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {*config.component_id},
                std::move(source_term_mesh_subset)));
        return ProcessLib::createNodalSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            source_term_mesh.getID(), variable_id, *config.component_id,
            parameters);
    }

    if (type == "Volumetric")
    {
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_source_term(
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {*config.component_id},
                std::move(source_term_mesh_subset)));
        return ProcessLib::createVolumetricSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            parameters, integration_order, shapefunction_order, variable_id,
            *config.component_id);
    }

    if (type == "Python")
    {
#ifdef OGS_USE_PYTHON
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table_source_term(
            dof_table_bulk.deriveBoundaryConstrainedMap(
                std::move(source_term_mesh_subset)));

        return ProcessLib::createPythonSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            source_term_mesh.getID(), variable_id, *config.component_id,
            integration_order, shapefunction_order,
            source_term_mesh.getDimension());
#else
        OGS_FATAL("OpenGeoSys has not been built with Python support.");
#endif
    }


    OGS_FATAL("Unknown source term type: `%s'.", type.c_str());
}
}  // namespace ProcessLib
