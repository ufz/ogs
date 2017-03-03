/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessOutput.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib
{

ProcessOutput::ProcessOutput(BaseLib::ConfigTree const& output_config)
{
    //! \ogs_file_param{prj__time_loop__processes__process__output__variables}
    auto const out_vars = output_config.getConfigSubtree("variables");

    //! \ogs_file_param{prj__time_loop__processes__process__output__variables__variable}
    for (auto out_var : out_vars.getConfigParameterList<std::string>("variable"))
    {
        if (output_variables.find(out_var) != output_variables.cend())
        {
            OGS_FATAL("output variable `%s' specified more than once.", out_var.c_str());
        }

        DBUG("adding output variable `%s'", out_var.c_str());
        output_variables.insert(out_var);
    }

    if (auto out_resid =
            //! \ogs_file_param{prj__time_loop__processes__process__output__output_extrapolation_residuals}
            output_config.getConfigParameterOptional<bool>("output_extrapolation_residuals"))
    {
        output_residuals = *out_resid;
    }
}


void doProcessOutput(
        std::string const& file_name,
        GlobalVector const& x,
        MeshLib::Mesh& mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<std::reference_wrapper<ProcessVariable>> const&
        process_variables,
        SecondaryVariableCollection secondary_variables,
        ProcessOutput const& process_output)
{
    DBUG("Process output.");

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

    // primary variables
    for (int variable_id = 0;
         variable_id < static_cast<int>(process_variables.size());
         ++variable_id)
    {
        ProcessVariable& pv = process_variables[variable_id];
        int const n_components = pv.getNumberOfComponents();
        global_component_offset = global_component_offset_next;
        global_component_offset_next += n_components;

        if (output_variables.find(pv.getName()) == output_variables.cend())
            continue;

        already_output.insert(pv.getName());

        DBUG("  process variable %s", pv.getName().c_str());

        auto& output_data = pv.getOrCreateMeshProperty();

        auto const num_comp = pv.getNumberOfComponents();

        for (int component_id = 0; component_id < num_comp; ++component_id)
        {
            auto const& mesh_subsets =
                dof_table.getMeshSubsets(variable_id,
                                                           component_id);
            for (auto const& mesh_subset : mesh_subsets)
            {
                auto const mesh_id = mesh_subset->getMeshID();
                for (auto const* node : mesh_subset->getNodes())
                {
                    MeshLib::Location const l(
                        mesh_id, MeshLib::MeshItemType::Node, node->getID());

                auto const global_component_id = global_component_offset + component_id;
                auto const index =
                        dof_table.getLocalIndex(
                            l, global_component_id, x.getRangeBegin(),
                            x.getRangeEnd());

                output_data[node->getID() * n_components + component_id] =
                        x_copy[index];
                }
            }
        }
    }

#ifndef USE_PETSC
    // the following section is for the output of secondary variables

    auto count_mesh_items = [](
        MeshLib::Mesh const& mesh, MeshLib::MeshItemType type) -> std::size_t
    {
        switch (type) {
        case MeshLib::MeshItemType::Cell: return mesh.getNumberOfElements();
        case MeshLib::MeshItemType::Node: return mesh.getNumberOfNodes();
        default: break; // avoid compiler warning
        }
        return 0;
    };

    auto get_or_create_mesh_property = [&mesh, &count_mesh_items](
        std::string const& property_name, MeshLib::MeshItemType type)
    {
        // Get or create a property vector for results.
        MeshLib::PropertyVector<double>* result = nullptr;

        auto const N = count_mesh_items(mesh, type);

        if (mesh.getProperties().existsPropertyVector<double>(property_name))
        {
            result = mesh.getProperties().template
                     getPropertyVector<double>(property_name);
        }
        else
        {
            result = mesh.getProperties().template
                     createNewPropertyVector<double>(property_name, type);
            result->resize(N);
        }
        assert(result && result->size() == N);

        return result;
    };

    auto add_secondary_var = [&](SecondaryVariable const& var,
                             std::string const& output_name)
    {
        assert(var.n_components == 1); // TODO implement other cases

        {
            DBUG("  secondary variable %s", output_name.c_str());

            auto result = get_or_create_mesh_property(
                              output_name, MeshLib::MeshItemType::Node);
            assert(result->size() == mesh.getNumberOfNodes());

            std::unique_ptr<GlobalVector> result_cache;
            auto const& nodal_values =
                    var.fcts.eval_field(x, dof_table, result_cache);

            // Copy result
            for (std::size_t i = 0; i < mesh.getNumberOfNodes(); ++i)
            {
                assert(!std::isnan(nodal_values[i]));
                (*result)[i] = nodal_values[i];
            }
        }

        if (process_output.output_residuals && var.fcts.eval_residuals)
        {
            DBUG("  secondary variable %s residual", output_name.c_str());
            auto const& property_name_res = output_name + "_residual";

            auto result = get_or_create_mesh_property(
                              property_name_res, MeshLib::MeshItemType::Cell);
            assert(result->size() == mesh.getNumberOfElements());

            std::unique_ptr<GlobalVector> result_cache;
            auto const& residuals =
                    var.fcts.eval_residuals(x, dof_table, result_cache);

            // Copy result
            for (std::size_t i = 0; i < mesh.getNumberOfElements(); ++i)
            {
                assert(!std::isnan(residuals[i]));
                (*result)[i] = residuals[i];
            }
        }
    };

    for (auto const& external_variable_name : output_variables)
    {
        if (!already_output.insert(external_variable_name).second) {
            // no insertion took place, output already done
            continue;
        }

        add_secondary_var(secondary_variables.get(external_variable_name),
                          external_variable_name);
    }

    // secondary variables output end
#else
    (void) secondary_variables;
#endif // USE_PETSC

    // Write output file
    DBUG("Writing output to \'%s\'.", file_name.c_str());
    MeshLib::IO::VtuInterface vtu_interface(&mesh, vtkXMLWriter::Binary, true);
    vtu_interface.writeToFile(file_name);
}

} // ProcessLib
