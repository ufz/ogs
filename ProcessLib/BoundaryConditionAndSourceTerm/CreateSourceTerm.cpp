/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSourceTerm.h"

#include "CreateNodalSourceTerm.h"
#include "CreateVolumetricSourceTerm.h"
#include "Python/CreatePythonSourceTerm.h"
#include "SourceTerm.h"
#include "SourceTermConfig.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table_bulk,
    const MeshLib::Mesh& source_term_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    [[maybe_unused]] std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process)
{
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    auto const type = config.config.peekConfigParameter<std::string>("type");

    // check basic data consistency
    if (variable_id >=
            static_cast<int>(dof_table_bulk.getNumberOfVariables()) ||
        config.component_id >=
            dof_table_bulk.getNumberOfVariableComponents(variable_id))
    {
        OGS_FATAL(
            "Variable id or component id too high. Actual values: ({:d}, "
            "{:d}), maximum values: ({:d}, {:d}).",
            variable_id, config.component_id,
            dof_table_bulk.getNumberOfVariables(),
            dof_table_bulk.getNumberOfVariableComponents(variable_id));
    }

    if (!source_term_mesh.getProperties()
             .template existsPropertyVector<std::size_t>(
                 MeshLib::getBulkIDString(MeshLib::MeshItemType::Node)))
    {
        OGS_FATAL(
            "The required bulk node ids map does not exist in the source term "
            "mesh '{:s}'.",
            source_term_mesh.getName());
    }
    std::vector<MeshLib::Node*> const& source_term_nodes =
        source_term_mesh.getNodes();
    DBUG(
        "Found {:d} nodes for source term at mesh '{:s}' for the variable {:d} "
        "and component {:d}",
        source_term_nodes.size(), source_term_mesh.getName(), variable_id,
        config.component_id);

    MeshLib::MeshSubset source_term_mesh_subset(source_term_mesh,
                                                source_term_nodes);

    if (type == "Nodal")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {config.component_id},
                std::move(source_term_mesh_subset));
        return ProcessLib::createNodalSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            source_term_mesh.getID(), variable_id, config.component_id,
            parameters);
    }

    if (type == "Line" || type == "Volumetric")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, {config.component_id},
                std::move(source_term_mesh_subset));
        auto const& bulk_mesh_dimension =
            dof_table_bulk.getMeshSubset(variable_id, config.component_id)
                .getMesh()
                .getDimension();
        return ProcessLib::createVolumetricSourceTerm(
            config.config, bulk_mesh_dimension, config.mesh,
            std::move(dof_table_source_term), parameters, integration_order,
            shapefunction_order);
    }

    if (type == "Python")
    {
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                std::move(source_term_mesh_subset));

        return ProcessLib::createPythonSourceTerm(
            config.config, config.mesh, std::move(dof_table_source_term),
            variable_id, config.component_id, integration_order,
            shapefunction_order, source_term_mesh.getDimension(),
            all_process_variables_for_this_process);
    }

    OGS_FATAL("Unknown source term type: `{:s}'.", type);
}
}  // namespace ProcessLib
