/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"

namespace ProcessLib
{
class Process;
class ProcessVariable;
class SecondaryVariableCollection;
struct IntegrationPointWriter;

/// Holds all data of a process that are needed for output.
class ProcessOutputData final
{
public:
    ProcessOutputData(
        std::vector<std::reference_wrapper<
            const std::vector<std::reference_wrapper<ProcessVariable>>>>&&
            process_variables_of_all_processes,
        const SecondaryVariableCollection& secondary_variables,
        const std::vector<std::unique_ptr<IntegrationPointWriter>>*
            integration_point_writers,
        std::vector<const NumLib::LocalToGlobalIndexMap*>&&
            bulk_mesh_dof_tables_of_all_processes,
        std::vector<const NumLib::LocalToGlobalIndexMap*>&&
            output_mesh_dof_tables_of_all_processes,
        std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>&&
            container_that_owns_output_mesh_dof_tables,
        MeshLib::Mesh& output_mesh);

    std::vector<std::reference_wrapper<ProcessVariable>> const&
    getProcessVariables(int const process_id) const
    {
        return process_variables_of_all_processes_[process_id].get();
    }

    SecondaryVariableCollection const& getSecondaryVariables() const
    {
        return secondary_variables_;
    }

    std::vector<std::unique_ptr<IntegrationPointWriter>> const*
    getIntegrationPointWriters() const
    {
        return integration_point_writers_;
    }

    NumLib::LocalToGlobalIndexMap const& getBulkMeshDofTable(
        int const process_id) const
    {
        return *bulk_mesh_dof_tables_of_all_processes_[process_id];
    }

    NumLib::LocalToGlobalIndexMap const& getOutputMeshDofTable(
        int const process_id) const
    {
        return *output_mesh_dof_tables_of_all_processes_[process_id];
    }

    std::vector<const NumLib::LocalToGlobalIndexMap*> const&
    getOutputMeshDofTablesOfAllProcesses() const
    {
        return output_mesh_dof_tables_of_all_processes_;
    }

    MeshLib::Mesh& getOutputMesh() const { return output_mesh_; }

private:
    /// Process variables of all processes.
    ///
    /// Each element of the container corresponds to a specific \c process_id of
    /// the Process.
    std::vector<std::reference_wrapper<
        const std::vector<std::reference_wrapper<ProcessVariable>>>>
        process_variables_of_all_processes_;

    SecondaryVariableCollection const& secondary_variables_;

    /// The list of integration point writers or \c nullptr if there are no
    /// integration point writers defined on the #output_mesh_.
    ///
    /// The latter is the case for output on submeshes.
    std::vector<std::unique_ptr<IntegrationPointWriter>> const*
        integration_point_writers_;

    /// D.o.f. tables for the full simulation domain of the Process this
    /// ProcessOutputData is associated with.
    ///
    /// Each element of the container corresponds to a specific \c process_id of
    /// the Process.
    std::vector<NumLib::LocalToGlobalIndexMap const*>
        bulk_mesh_dof_tables_of_all_processes_;

    /// D.o.f tables for the given #output_mesh_.
    ///
    /// In the case of submesh output these d.o.f. tables are different from the
    /// #bulk_mesh_dof_tables_of_all_processes_.
    ///
    /// Each element of the container corresponds to a specific \c process_id of
    /// the Process this ProcessOutputData is associated with.
    std::vector<NumLib::LocalToGlobalIndexMap const*>
        output_mesh_dof_tables_of_all_processes_;

    /// Actual data backing the pointers in
    /// #output_mesh_dof_tables_of_all_processes_.
    ///
    /// This container is populated in the case of submesh output only.
    std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>
        container_that_owns_output_mesh_dof_tables_;

    /// The mesh to which output shall be written.
    ///
    /// This can be the entire simulation domain of the Process this
    /// ProcessOutputData is associated with or a submesh thereof.
    MeshLib::Mesh& output_mesh_;
};

/// Extracts data necessary for output from the given \c process.
ProcessOutputData createProcessOutputData(Process const& process,
                                          std::size_t const n_processes,
                                          MeshLib::Mesh& output_mesh);

}  // namespace ProcessLib
