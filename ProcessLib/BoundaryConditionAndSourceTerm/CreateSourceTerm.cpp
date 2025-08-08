/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateSourceTerm.h"

#include <cassert>
#include <numeric>

#include "CreateAnchorTerm.h"
#include "CreateEmbeddedAnchor.h"
#include "CreateNodalSourceTerm.h"
#include "CreateVolumetricSourceTerm.h"
#include "Python/CreatePythonSourceTerm.h"
#include "SourceTerm.h"
#include "SourceTermConfig.h"

namespace ProcessLib
{
std::unique_ptr<SourceTermBase> createSourceTerm(
    const SourceTermConfig& config,
    const NumLib::LocalToGlobalIndexMap& dof_table_bulk,
    const MeshLib::Mesh& source_term_mesh, const int variable_id,
    const unsigned integration_order, const unsigned shapefunction_order,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    [[maybe_unused]] std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process,
    const MeshLib::Mesh& bulk_mesh)
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
                 MeshLib::getBulkIDString(MeshLib::MeshItemType::Node)) &&
        (type != "EmbeddedAnchor"))
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
    if (type == "Anchor")
    {
        const int number_of_components =
            dof_table_bulk.getNumberOfVariableComponents(variable_id);
        std::vector<int> component_ids(number_of_components);
        std::iota(std::begin(component_ids), std::end(component_ids), 0);
        auto dof_table_source_term =
            dof_table_bulk.deriveBoundaryConstrainedMap(
                variable_id, component_ids, std::move(source_term_mesh_subset));
        const int bulk_mesh_dimension =
            dof_table_bulk.getMeshSubset(variable_id, 0)
                .getMesh()
                .getDimension();
        if (bulk_mesh_dimension != number_of_components)
        {
            OGS_FATAL(
                "For the Anchor source term type,"
                "the bulk mesh dimension needs to be the same "
                "as the number of process variable components.");
        }
        switch (bulk_mesh_dimension)
        {
            case 2:
                return ProcessLib::createAnchorTerm<2>(
                    config.config, config.mesh,
                    std::move(dof_table_source_term), source_term_mesh.getID(),
                    variable_id, parameters);
            case 3:
                return ProcessLib::createAnchorTerm<3>(
                    config.config, config.mesh,
                    std::move(dof_table_source_term), source_term_mesh.getID(),
                    variable_id, parameters);
            default:
                OGS_FATAL(
                    "Anchor can not be instantiated "
                    "for mesh dimensions other than two or three. "
                    "{}-dimensional mesh was given.",
                    bulk_mesh_dimension);
        }
    }
    if (type == "EmbeddedAnchor")
    {
        const int number_of_components =
            dof_table_bulk.getNumberOfVariableComponents(variable_id);

        const int bulk_mesh_dimension = bulk_mesh.getDimension();
        if (bulk_mesh_dimension != number_of_components)
        {
            OGS_FATAL(
                "For the EmbeddedAnchor source term type,"
                "the bulk mesh dimension needs to be the same "
                "as the number of process variable components.");
        }
        switch (bulk_mesh_dimension)
        {
            case 2:
                return ProcessLib::createEmbeddedAnchor<2>(
                    config.config, config.mesh, bulk_mesh, dof_table_bulk,
                    source_term_mesh.getID(), variable_id);
            case 3:
                return ProcessLib::createEmbeddedAnchor<3>(
                    config.config, config.mesh, bulk_mesh, dof_table_bulk,
                    source_term_mesh.getID(), variable_id);
            default:
                OGS_FATAL(
                    "Anchor can not be instantiated "
                    "for mesh dimensions other than two or three. "
                    "{}-dimensional mesh was given.",
                    bulk_mesh_dimension);
        }
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
