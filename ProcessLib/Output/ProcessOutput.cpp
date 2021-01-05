/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessOutput.h"

#ifndef _WIN32
#ifndef __APPLE__
#include <cfenv>
#endif  // __APPLE__
#endif  // _WIN32

#include "InfoLib/GitInfo.h"
#include "IntegrationPointWriter.h"
#include "MathLib/LinAlg/LinAlg.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

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
#ifdef USE_PETSC
    std::size_t const global_vector_size =
        nodal_values.getLocalSize() + nodal_values.getGhostSize();
#else
    std::size_t const global_vector_size = nodal_values.size();
#endif
    if (nodal_values_mesh.size() != global_vector_size)
    {
        OGS_FATAL(
            "Secondary variable `{:s}' did not evaluate to the right "
            "number of components. Expected: {:d}, actual: {:d}.",
            var.name, nodal_values_mesh.size(), global_vector_size);
    }

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
            "components. "
            "Expected: {:d}, actual: {:d}",
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

namespace ProcessLib
{
void addProcessDataToMesh(
    const double t, std::vector<GlobalVector*> const& x, int const process_id,
    MeshLib::Mesh& mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
    SecondaryVariableCollection const& secondary_variables,
    bool const output_secondary_variable,
    std::vector<std::unique_ptr<IntegrationPointWriter>> const* const
        integration_point_writer,
    OutputDataSpecification const& process_output)
{
    DBUG("Process output data.");

    addOgsVersion(mesh);

    // Copy result
#ifdef USE_PETSC
    // TODO It is also possible directly to copy the data for single process
    // variable to a mesh property. It needs a vector of global indices and
    // some PETSc magic to do so.
    std::vector<double> x_copy(x[process_id]->getLocalSize() +
                               x[process_id]->getGhostSize());
#else
    std::vector<double> x_copy(x[process_id]->size());
#endif
    x[process_id]->copyValues(x_copy);

    auto const& output_variables = process_output.output_variables;
    std::set<std::string> already_output;

    int global_component_offset = 0;
    int global_component_offset_next = 0;

    const auto number_of_dof_variables =
        dof_table[process_id]->getNumberOfVariables();
    // primary variables
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        int const n_components = pv.getNumberOfGlobalComponents();
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

        DBUG("  process variable {:s}", pv.getName());

        auto const num_comp = pv.getNumberOfGlobalComponents();
        auto& output_data = *MeshLib::getOrCreateMeshProperty<double>(
            mesh, pv.getName(), MeshLib::MeshItemType::Node, num_comp);

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subset = dof_table[process_id]->getMeshSubset(
                sub_meshset_id, component_id);
            auto const mesh_id = mesh_subset.getMeshID();
            for (auto const* node : mesh_subset.getNodes())
            {
                MeshLib::Location const l(mesh_id, MeshLib::MeshItemType::Node,
                                          node->getID());

                auto const global_component_id =
                    global_component_offset + component_id;
                auto const index = dof_table[process_id]->getLocalIndex(
                    l, global_component_id, x[process_id]->getRangeBegin(),
                    x[process_id]->getRangeEnd());

                output_data[node->getID() * n_components + component_id] =
                    x_copy[index];
            }
        }
    }

    if (output_secondary_variable)
    {
        for (auto const& external_variable_name : secondary_variables)
        {
            auto const& name = external_variable_name.first;
            if (!already_output.insert(name).second)
            {
                // no insertion took place, output already done
                continue;
            }

            addSecondaryVariableNodes(
                t, x, dof_table, secondary_variables.get(name), name, mesh);

            if (process_output.output_residuals)
            {
                addSecondaryVariableResiduals(
                    t, x, dof_table, secondary_variables.get(name), name, mesh);
            }
        }
    }

    if (integration_point_writer != nullptr)
    {
        addIntegrationPointWriter(mesh, *integration_point_writer);
    }
}

void makeOutput(std::string const& file_name, MeshLib::Mesh const& mesh,
                bool const compress_output, int const data_mode,
                OutputType const file_type, int const time_step,
                double const time)
{
    // Write output file
    DBUG("Writing output to '{:s}'.", file_name);

    // Store floating-point exception handling. Output of NaN's triggers
    // floating point exceptions. Because we are not debugging VTK (or other
    // libraries) at this point, the possibly set exceptions are temporary
    // disabled and restored before end of the function.
#ifndef _WIN32
#ifndef __APPLE__
    fenv_t fe_env;
    fegetenv(&fe_env);
    fesetenv(FE_DFL_ENV);  // Set default environment effectively disabling
                           // exceptions.
#endif  //_WIN32
#endif  //__APPLE__

    switch (file_type)
    {
        case OutputType::vtk:
        {
            MeshLib::IO::VtuInterface vtu_interface(&mesh, data_mode,
                                                    compress_output);
            vtu_interface.writeToFile(file_name);
            break;
        }

        case OutputType::xdmf:
        {
#ifdef OGS_USE_XDMF
            MeshLib::IO::XdmfHdfWriter writer(
                mesh, file_name, time_step);
            writer.writeStep(time_step, time);
#else
            // silence compiler warnings
            (void)(time_step);
            (void)(time);
            OGS_FATAL(
                "Trying to write Xdmf file but OGS was not built with "
                "Xdmf-support.");
#endif
        }
    }

        // Restore floating-point exception handling.
#ifndef _WIN32
#ifndef __APPLE__
    fesetenv(&fe_env);
#endif  //_WIN32
#endif  //__APPLE__
}

}  // namespace ProcessLib
