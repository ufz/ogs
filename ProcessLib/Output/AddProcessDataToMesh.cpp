/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AddProcessDataToMesh.h"

#include "InfoLib/GitInfo.h"
#include "IntegrationPointWriter.h"
#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Process.h"

/// Copies the ogs_version string containing the release number and the git
/// hash.
static void addOgsVersion(MeshLib::Mesh& mesh)
{
    auto& ogs_version_field = *MeshLib::getOrCreateMeshProperty<char>(
        mesh, GitInfoLib::GitInfo::OGS_VERSION,
        MeshLib::MeshItemType::IntegrationPoint, 1);

    ogs_version_field.assign(GitInfoLib::GitInfo::ogs_version.begin(),
                             GitInfoLib::GitInfo::ogs_version.end());
}

static void addSecondaryVariableNodes(
    double const t,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    ProcessLib::SecondaryVariable const& var,
    std::string const& output_name,
    MeshLib::Mesh& mesh)
{
    DBUG("  secondary variable {:s}", output_name);

    auto& nodal_values_mesh = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, output_name, MeshLib::MeshItemType::Node,
        var.fcts.num_components);
    if (nodal_values_mesh.size() !=
        mesh.getNumberOfNodes() * var.fcts.num_components)
    {
        OGS_FATAL(
            "Nodal property `{:s}' does not have the right number of "
            "components. Expected: {:d}, actual: {:d}",
            output_name,
            mesh.getNumberOfNodes() * var.fcts.num_components,
            nodal_values_mesh.size());
    }

    std::unique_ptr<GlobalVector> result_cache;
    auto const& nodal_values =
        var.fcts.eval_field(t, x, dof_table, result_cache);

    // Copy result
    nodal_values.copyValues(nodal_values_mesh);
}

static void addSecondaryVariableResiduals(
    double const t,
    std::vector<GlobalVector*> const& x,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    ProcessLib::SecondaryVariable const& var,
    std::string const& output_name,
    MeshLib::Mesh& mesh)
{
    if (!var.fcts.eval_residuals)
    {
        return;
    }

    DBUG("  secondary variable {:s} residual", output_name);
    auto const& property_name_res = output_name + "_residual";

    auto& residuals_mesh = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, property_name_res, MeshLib::MeshItemType::Cell,
        var.fcts.num_components);
    if (residuals_mesh.size() !=
        mesh.getNumberOfElements() * var.fcts.num_components)
    {
        OGS_FATAL(
            "Cell property `{:s}' does not have the right number of "
            "components. Expected: {:d}, actual: {:d}",
            property_name_res,
            mesh.getNumberOfElements() * var.fcts.num_components,
            residuals_mesh.size());
    }

    std::unique_ptr<GlobalVector> result_cache;
    auto const& residuals =
        var.fcts.eval_residuals(t, x, dof_table, result_cache);
#ifdef USE_PETSC
    std::size_t const global_vector_size =
        residuals.getLocalSize() + residuals.getGhostSize();
#else
    std::size_t const global_vector_size = residuals.size();
#endif
    if (residuals_mesh.size() != global_vector_size)
    {
        OGS_FATAL(
            "The residual of secondary variable `{:s}' did not evaluate to the "
            "right number of components. Expected: {:d}, actual: {:d}.",
            var.name, residuals_mesh.size(), global_vector_size);
    }

    // Copy result
    residuals.copyValues(residuals_mesh);
}

static std::vector<double> copySolutionVector(GlobalVector const& x)
{
    std::vector<double> x_copy;
    x.copyValues(x_copy);
    return x_copy;
}

MeshLib::PropertyVector<std::size_t> const* getBulkNodeIdMapForPetscIfNecessary(
    [[maybe_unused]] MeshLib::Mesh const& mesh,
    [[maybe_unused]] NumLib::LocalToGlobalIndexMap const& mesh_dof_table,
    [[maybe_unused]] NumLib::LocalToGlobalIndexMap const& full_mesh_dof_table)
{
#ifdef USE_PETSC

    if (&full_mesh_dof_table != &mesh_dof_table)
    {
        if (!mesh.getProperties().existsPropertyVector<std::size_t>(
                "bulk_node_ids"))
        {
            OGS_FATAL(
                "The required bulk node ids map does not exist in "
                "the boundary mesh '{:s}' or has the wrong data "
                "type (should be equivalent to C++ data type "
                "std::size_t which is an unsigned integer of size "
                "{:d} or UInt64 in vtk terminology).",
                mesh.getName(), sizeof(std::size_t));
        }
        return mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids");
    }
#endif

    return nullptr;
}

static GlobalIndexType getIndexForComponentInSolutionVector(
    std::size_t const mesh_id, std::size_t const node_id,
    [[maybe_unused]] bool const is_ghost_node, int const global_component_id,
    GlobalVector const& x, NumLib::LocalToGlobalIndexMap const& mesh_dof_table,
    [[maybe_unused]] NumLib::LocalToGlobalIndexMap const& full_mesh_dof_table,
    [[maybe_unused]] MeshLib::PropertyVector<std::size_t> const* const
        bulk_node_id_map)
{
#ifdef USE_PETSC
    if (is_ghost_node && &full_mesh_dof_table != &mesh_dof_table)
    {
        auto const bulk_node_id = (*bulk_node_id_map)[node_id];
        std::size_t const bulk_mesh_id = 0;

        MeshLib::Location const loc(bulk_mesh_id, MeshLib::MeshItemType::Node,
                                    bulk_node_id);

        // early return!
        return full_mesh_dof_table.getLocalIndex(
            loc, global_component_id, x.getRangeBegin(), x.getRangeEnd());
    }
#endif

    MeshLib::Location const loc(mesh_id, MeshLib::MeshItemType::Node, node_id);

    return mesh_dof_table.getLocalIndex(loc, global_component_id,
                                        x.getRangeBegin(), x.getRangeEnd());
}

static bool isGhostNode([[maybe_unused]] MeshLib::Mesh const& mesh,
                        [[maybe_unused]] std::size_t node_id)
{
#ifndef USE_PETSC
    return false;
#else
    return static_cast<MeshLib::NodePartitionedMesh const&>(mesh).isGhostNode(
        node_id);
#endif
}

/**
 * Adds data for all given primary variables to the given \c mesh.
 *
 * \param mesh the mesh the data is added to.
 * \param x the global solution vector providing the data.
 * \param process_variables the primary variables comprising \c x.
 * \param output_variables the names of variables that can be added to the
 *        \c mesh.
 * \param mesh_dof_table the d.o.f. table related to the passed \c mesh and
 *        solution vector \c x.
 * \param full_mesh_dof_table the d.o.f. table related to the full simulation
 *        domain and the passed solution vector \c x.
 *
 * \note Usually \c mesh and the full simulation domain are the same. But if
 *       output should be written to a sub mesh, they will differ. In that case,
 *       also the two d.o.f. table will be different from each other.
 *
 * \return The names of all variables that have been written to the \c mesh.
 */
static std::set<std::string> addPrimaryVariablesToMesh(
    MeshLib::Mesh& mesh,
    GlobalVector const& x,
    std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>> const&
        process_variables,
    std::set<std::string> const& output_variables,
    NumLib::LocalToGlobalIndexMap const& mesh_dof_table,
    NumLib::LocalToGlobalIndexMap const& full_mesh_dof_table)
{
    auto const x_copy = copySolutionVector(x);
    std::set<std::string> names_of_already_output_variables;

    const auto number_of_dof_variables = mesh_dof_table.getNumberOfVariables();
    assert(number_of_dof_variables == process_variables.size());

    int global_component_offset = 0;
    int global_component_offset_next = 0;

    auto const* const bulk_node_id_map = getBulkNodeIdMapForPetscIfNecessary(
        mesh, mesh_dof_table, full_mesh_dof_table);

    for (int variable_id = 0; variable_id < number_of_dof_variables;
         ++variable_id)
    {
        auto const& pv = process_variables[variable_id].get();
        auto const n_components = pv.getNumberOfGlobalComponents();

        // increase global component offset even if we do not add anything in
        // this iteration
        global_component_offset = global_component_offset_next;
        global_component_offset_next += n_components;

        if (output_variables.find(pv.getName()) == output_variables.cend())
        {
            continue;
        }

        names_of_already_output_variables.insert(pv.getName());

        DBUG("  process variable {:s}", pv.getName());

        auto& output_data = *MeshLib::getOrCreateMeshProperty<double>(
            mesh, pv.getName(), MeshLib::MeshItemType::Node, n_components);

        // mesh subsets are the same for all components
        int const dummy_component_id = 0;
        auto const& mesh_subset =
            mesh_dof_table.getMeshSubset(variable_id, dummy_component_id);
        auto const mesh_id = mesh_subset.getMeshID();

        for (auto const* node : mesh_subset.getNodes())
        {
            auto const node_id = node->getID();
            auto const is_ghost_node = isGhostNode(mesh, node_id);

            for (int component_id = 0; component_id < n_components;
                 ++component_id)
            {
                auto const global_component_id =
                    global_component_offset + component_id;

                auto const in_index = getIndexForComponentInSolutionVector(
                    mesh_id, node_id, is_ghost_node, global_component_id, x,
                    mesh_dof_table, full_mesh_dof_table, bulk_node_id_map);

                // per node ordering of components
                auto const out_index = node_id * n_components + component_id;

                output_data[out_index] = x_copy[in_index];
            }
        }
    }

    return names_of_already_output_variables;
}

static void addSecondaryVariablesToMesh(
    ProcessLib::SecondaryVariableCollection const& secondary_variables,
    std::set<std::string>& names_of_already_output_variables, const double t,
    std::vector<GlobalVector*> const& xs, MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_tables,
    bool const output_residuals)
{
    for (auto const& external_variable_name : secondary_variables)
    {
        auto const& name = external_variable_name.first;
        if (!names_of_already_output_variables.insert(name).second)
        {
            // no insertion took place, output already done
            continue;
        }

        addSecondaryVariableNodes(t, xs, dof_tables,
                                  secondary_variables.get(name), name, mesh);

        if (output_residuals)
        {
            addSecondaryVariableResiduals(
                t, xs, dof_tables, secondary_variables.get(name), name, mesh);
        }
    }
}

namespace ProcessLib
{
void addProcessDataToMesh(
    const double t, std::vector<GlobalVector*> const& xs, int const process_id,
    MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& mesh_dof_tables,
    Process const& process, bool const output_secondary_variables,
    OutputDataSpecification const& process_output)
{
    DBUG("Process output data.");

    auto const& process_variables = process.getProcessVariables(process_id);
    auto const& secondary_variables = process.getSecondaryVariables();
    auto const* const integration_point_writers =
        process.getIntegrationPointWriter(mesh);
    auto const& output_variables = process_output.output_variables;
    auto const& full_mesh_dof_table = process.getDOFTable(process_id);

    addOgsVersion(mesh);

    auto names_of_already_output_variables = addPrimaryVariablesToMesh(
        mesh, *xs[process_id], process_variables, output_variables,
        *mesh_dof_tables[process_id], full_mesh_dof_table);

    if (output_secondary_variables)
    {
        addSecondaryVariablesToMesh(
            secondary_variables, names_of_already_output_variables, t, xs, mesh,
            mesh_dof_tables, process_output.output_residuals);
    }

    if (integration_point_writers)
    {
        addIntegrationPointDataToMesh(mesh, *integration_point_writers);
    }
}

void addProcessDataToMesh(const double t, std::vector<GlobalVector*> const& xs,
                          int const process_id, MeshLib::Mesh& mesh,
                          Process const& process,
                          bool const output_secondary_variables,
                          OutputDataSpecification const& process_output)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> full_mesh_dof_tables(
        xs.size());
    for (std::size_t i = 0; i < xs.size(); ++i)
    {
        full_mesh_dof_tables[i] = &process.getDOFTable(i);
    }

    addProcessDataToMesh(t, xs, process_id, mesh, full_mesh_dof_tables, process,
                         output_secondary_variables, process_output);
}
}  // namespace ProcessLib
