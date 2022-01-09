/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ProcessOutputData.h"

#include "ProcessLib/Process.h"

namespace
{
/// Checks if the given \c mesh is the simulation domain of the given \c
/// process.
bool isSimulationDomain(MeshLib::Mesh const& mesh,
                        ProcessLib::Process const& process)
{
    return mesh == process.getMesh();
}

std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>
computeDofTablesForSubmesh(ProcessLib::Process const& process,
                           MeshLib::Mesh const& submesh,
                           std::size_t const n_processes)
{
    std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>
        submesh_dof_tables;
    submesh_dof_tables.reserve(n_processes);

    for (std::size_t i = 0; i < n_processes; ++i)
    {
        submesh_dof_tables.push_back(
            process.getDOFTable(i).deriveBoundaryConstrainedMap(
                MeshLib::MeshSubset{submesh, submesh.getNodes()}));
    }

    return submesh_dof_tables;
}

std::vector<NumLib::LocalToGlobalIndexMap const*> toNonOwning(
    std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>> const&
        dof_tables)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_table_pointers;

    dof_table_pointers.reserve(dof_tables.size());
    transform(cbegin(dof_tables), cend(dof_tables),
              back_inserter(dof_table_pointers),
              [](std::unique_ptr<NumLib::LocalToGlobalIndexMap> const& p)
              { return p.get(); });

    return dof_table_pointers;
}

std::vector<NumLib::LocalToGlobalIndexMap const*> getDofTablesOfAllProcesses(
    ProcessLib::Process const& process, std::size_t const n_processes)
{
    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables_of_all_procs(
        n_processes);

    for (std::size_t proc_id = 0; proc_id < n_processes; ++proc_id)
    {
        dof_tables_of_all_procs[proc_id] = &process.getDOFTable(proc_id);
    }

    return dof_tables_of_all_procs;
}

/// Computes the d.o.f. tables for the given \c output_mesh.
///
/// These are the passed \c bulk_mesh_dof_tables for output of the entire
/// simulation domain of the given \c process. In the case of submesh output
/// d.o.f. tables for the submesh will be computed.
///
/// \return A pair of (vector of d.o.f. table pointers, vector of d.o.f. table
/// storage), where the latter is populated in the case of submesh output only.
///
/// Each element in the returned vectors corresponds to a specific \c process_id
/// of the \c process.
decltype(auto) computeOutputMeshDofTables(
    ProcessLib::Process const& process,
    MeshLib::Mesh const& output_mesh,
    std::vector<NumLib::LocalToGlobalIndexMap const*> const&
        bulk_mesh_dof_tables)
{
    if (isSimulationDomain(output_mesh, process))
    {
        return std::pair(
            bulk_mesh_dof_tables /* will be copied */,
            std::vector<std::unique_ptr<NumLib::LocalToGlobalIndexMap>>{});
    }

    auto const n_processes = bulk_mesh_dof_tables.size();

    // TODO Currently these d.o.f. tables will be recomputed everytime we write
    // output. That should be avoided in the future.
    auto container_that_owns_output_mesh_dof_tables =
        computeDofTablesForSubmesh(process, output_mesh, n_processes);

    auto output_mesh_dof_tables =
        toNonOwning(container_that_owns_output_mesh_dof_tables);

    return std::pair(std::move(output_mesh_dof_tables),
                     std::move(container_that_owns_output_mesh_dof_tables));
}

std::vector<std::reference_wrapper<
    const std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>>>
getProcessVariablesOfAllProcesses(ProcessLib::Process const& process,
                                  std::size_t const n_processes)
{
    std::vector<std::reference_wrapper<
        const std::vector<std::reference_wrapper<ProcessLib::ProcessVariable>>>>
        pvs_of_all_procs;
    pvs_of_all_procs.reserve(n_processes);

    for (std::size_t proc_id = 0; proc_id < n_processes; ++proc_id)
    {
        pvs_of_all_procs.emplace_back(process.getProcessVariables(proc_id));
    }

    return pvs_of_all_procs;
}

std::vector<std::unique_ptr<ProcessLib::IntegrationPointWriter>> const*
getIntegrationPointWriters(ProcessLib::Process const& process,
                           MeshLib::Mesh const& output_mesh)
{
    return isSimulationDomain(output_mesh, process)
               ? &process.getIntegrationPointWriters()
               : nullptr;
}

}  // namespace

namespace ProcessLib
{

ProcessOutputData createProcessOutputData(Process const& process,
                                          std::size_t const n_processes,
                                          MeshLib::Mesh& output_mesh)
{
    auto bulk_mesh_dof_tables =
        getDofTablesOfAllProcesses(process, n_processes);

    auto [output_mesh_dof_tables, container_that_owns_output_mesh_dof_tables] =
        computeOutputMeshDofTables(process, output_mesh, bulk_mesh_dof_tables);

    return {getProcessVariablesOfAllProcesses(process, n_processes),
            process.getSecondaryVariables(),
            ::getIntegrationPointWriters(process, output_mesh),
            std::move(bulk_mesh_dof_tables),
            std::move(output_mesh_dof_tables),
            std::move(container_that_owns_output_mesh_dof_tables),
            output_mesh};
}

ProcessOutputData::ProcessOutputData(
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
    MeshLib::Mesh& output_mesh)
    : process_variables_of_all_processes_(
          std::move(process_variables_of_all_processes)),
      secondary_variables_(secondary_variables),
      integration_point_writers_(integration_point_writers),
      bulk_mesh_dof_tables_of_all_processes_(
          std::move(bulk_mesh_dof_tables_of_all_processes)),
      output_mesh_dof_tables_of_all_processes_(
          std::move(output_mesh_dof_tables_of_all_processes)),
      container_that_owns_output_mesh_dof_tables_(
          std::move(container_that_owns_output_mesh_dof_tables)),
      output_mesh_(output_mesh)
{
    auto const n_proc_pvs = process_variables_of_all_processes_.size();
    auto const n_proc_bulk = bulk_mesh_dof_tables_of_all_processes_.size();
    auto const n_proc_out = output_mesh_dof_tables_of_all_processes_.size();
    auto const n_proc_own = container_that_owns_output_mesh_dof_tables_.size();

    if (n_proc_pvs != n_proc_bulk)
    {
        OGS_FATAL(
            "Mismatch in number of processes (PVs vs. bulk mesh d.o.f. "
            "tables): {} != {}",
            n_proc_pvs, n_proc_bulk);
    }

    if (n_proc_pvs != n_proc_out)
    {
        OGS_FATAL(
            "Mismatch in number of processes (PVs vs. output mesh d.o.f. "
            "tables): {} != {}",
            n_proc_pvs, n_proc_out);
    }

    // n_proc_own is nonzero only for submesh output
    if (n_proc_own != 0 && n_proc_pvs != n_proc_own)
    {
        OGS_FATAL(
            "Mismatch in number of processes (PVs vs. output mesh d.o.f. "
            "tables, owning): {} != {}",
            n_proc_pvs, n_proc_own);
    }
}

}  // namespace ProcessLib
