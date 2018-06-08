/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessOutput.h"

#include "BaseLib/BuildInfo.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

#include "IntegrationPointWriter.h"


/// Copies the ogs_version string containing the release number and the git
/// hash.
static void addOgsVersion(MeshLib::Mesh& mesh)
{
    auto& ogs_version_field = *MeshLib::getOrCreateMeshProperty<char>(
        mesh, "OGS_VERSION", MeshLib::MeshItemType::IntegrationPoint, 1);

    ogs_version_field.assign(BaseLib::BuildInfo::ogs_version.begin(),
                             BaseLib::BuildInfo::ogs_version.end());
}

static void addSecondaryVariableNodes(
    double const t,
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    ProcessLib::SecondaryVariable const& var,
    std::string const& output_name,
    MeshLib::Mesh& mesh)
{
    DBUG("  secondary variable %s", output_name.c_str());

    auto& nodal_values_mesh = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, output_name, MeshLib::MeshItemType::Node,
        var.fcts.num_components);
    if (nodal_values_mesh.size() !=
        mesh.getNumberOfNodes() * var.fcts.num_components)
    {
        OGS_FATAL(
            "Nodal property `%s' does not have the right number of "
            "components. Expected: %d, actual: %d",
            output_name.c_str(),
            mesh.getNumberOfNodes() * var.fcts.num_components,
            nodal_values_mesh.size());
    }

    std::unique_ptr<GlobalVector> result_cache;
    auto const& nodal_values =
        var.fcts.eval_field(t, x, dof_table, result_cache);
#ifdef USE_PETSC
    std::size_t const global_vector_size =
        nodal_values.getLocalSize() + nodal_values.getGhostSize();
#else
    std::size_t const global_vector_size = nodal_values.size();
#endif
    if (nodal_values_mesh.size() != global_vector_size)
    {
        OGS_FATAL(
            "Secondary variable `%s' did not evaluate to the right "
            "number of components. Expected: %d, actual: %d.",
            var.name.c_str(), nodal_values_mesh.size(), global_vector_size);
    }

    // Copy result
    nodal_values.copyValues(nodal_values_mesh);
}

static void addSecondaryVariableResiduals(
    double const t,
    GlobalVector const& x,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    ProcessLib::SecondaryVariable const& var,
    std::string const& output_name,
    MeshLib::Mesh& mesh)
{
    if (!var.fcts.eval_residuals)
        return;

    DBUG("  secondary variable %s residual", output_name.c_str());
    auto const& property_name_res = output_name + "_residual";

    auto& residuals_mesh = *MeshLib::getOrCreateMeshProperty<double>(
        mesh, property_name_res, MeshLib::MeshItemType::Cell,
        var.fcts.num_components);
    if (residuals_mesh.size() !=
        mesh.getNumberOfElements() * var.fcts.num_components)
    {
        OGS_FATAL(
            "Cell property `%s' does not have the right number of components. "
            "Expected: %d, actual: %d",
            property_name_res.c_str(),
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
            "The residual of secondary variable `%s' did not evaluate to the "
            "right number of components. Expected: %d, actual: %d.",
            var.name.c_str(), residuals_mesh.size(), global_vector_size);
    }

    // Copy result
    residuals.copyValues(residuals_mesh);
}

namespace ProcessLib
{
void processOutputData(
    const double t,
    GlobalVector const& x,
    MeshLib::Mesh& mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    SecondaryVariableCollection secondary_variables,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
        integration_point_writer,
    ProcessOutput const& process_output)
{
    DBUG("Process output data.");

    addOgsVersion(mesh);

    // Copy result
#ifdef USE_PETSC
    // TODO It is also possible directly to copy the data for single process
    // variable to a mesh property. It needs a vector of global indices and
    // some PETSc magic to do so.
    std::vector<double> x_copy(x.getLocalSize() + x.getGhostSize());
#else
    std::vector<double> x_copy(x.size());
#endif
    x.copyValues(x_copy);

    auto const& output_variables = process_output.output_variables;
    std::set<std::string> already_output;

    int global_component_offset = 0;
    int global_component_offset_next = 0;

    const auto number_of_dof_variables = dof_table.getNumberOfVariables();
    // primary variables
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        int const n_components = pv.getNumberOfComponents();
        // If (number_of_dof_variables==1), the case is either the staggered
        // scheme being applied or a single PDE being solved.
        const int sub_meshset_id =
            (number_of_dof_variables == 1) ? 0 : variable_id;

        if (number_of_dof_variables > 1)
        {
            global_component_offset = global_component_offset_next;
            global_component_offset_next += n_components;
        }

        if (output_variables.find(pv.getName()) == output_variables.cend())
        {
            continue;
        }

        already_output.insert(pv.getName());

        DBUG("  process variable %s", pv.getName().c_str());

        auto& output_data = pv.getOrCreateMeshProperty();

        auto const num_comp = pv.getNumberOfComponents();

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subset =
                dof_table.getMeshSubset(sub_meshset_id, component_id);
            auto const mesh_id = mesh_subset.getMeshID();
            for (auto const* node : mesh_subset.getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                          node->getID());

                auto const global_component_id =
                    global_component_offset + component_id;
                auto const index = dof_table.getLocalIndex(
                    l, global_component_id, x.getRangeBegin(), x.getRangeEnd());

                output_data[node->getID() * n_components + component_id] =
                    x_copy[index];
            }
        }
    }

    // Secondary variables output
    for (auto const& external_variable_name : output_variables)
    {
        if (!already_output.insert(external_variable_name).second)
        {
            // no insertion took place, output already done
            continue;
        }

        addSecondaryVariableNodes(
            t, x, dof_table, secondary_variables.get(external_variable_name),
            external_variable_name, mesh);
        if (process_output.output_residuals)
        {
            addSecondaryVariableResiduals(
                t, x, dof_table,
                secondary_variables.get(external_variable_name),
                external_variable_name, mesh);
        }
    }

    addIntegrationPointWriter(mesh, integration_point_writer);
}

void makeOutput(std::string const& file_name, MeshLib::Mesh& mesh,
                bool const compress_output, int const data_mode)
{
    // Write output file
    DBUG("Writing output to \'%s\'.", file_name.c_str());
    MeshLib::IO::VtuInterface vtu_interface(&mesh, data_mode, compress_output);
    vtu_interface.writeToFile(file_name);
}

}  // namespace ProcessLib
