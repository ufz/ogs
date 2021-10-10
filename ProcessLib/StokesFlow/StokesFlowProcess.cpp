/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StokesFlowProcess.h"

#include <cassert>

#include "CreateLocalAssemblers.h"
#include "MeshLib/Elements/Utils.h"

namespace ProcessLib
{
namespace StokesFlow
{
template <int GlobalDim>
StokesFlowProcess<GlobalDim>::StokesFlowProcess(
    std::string name,
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    StokesFlowProcessData&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    bool const use_monolithic_scheme)
    : Process(std::move(name), mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), use_monolithic_scheme),
      _process_data(std::move(process_data))
{
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, _base_nodes);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubset> all_mesh_subsets_single_component{
        *_mesh_subset_all_nodes};
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    assert(_use_monolithic_scheme);
    {
        // For vector variables, in this case liquid velocity.
        int const process_id = 0;
        std::vector<MeshLib::MeshSubset> all_mesh_subsets(
            getProcessVariables(process_id)[0]
                .get()
                .getNumberOfGlobalComponents(),
            *_mesh_subset_all_nodes);

        // For scalar variables, in this case pressure.
        all_mesh_subsets.push_back(*_mesh_subset_base_nodes);

        std::vector<int> const vec_n_components{GlobalDim, 1};

        _local_to_global_index_map =
            std::make_unique<NumLib::LocalToGlobalIndexMap>(
                std::move(all_mesh_subsets), vec_n_components,
                NumLib::ComponentOrder::BY_LOCATION);
        assert(_local_to_global_index_map);
    }
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    int const process_id = 0;
    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Todo: may move LocalDataInitializer.h and CreateLocalAssemblers.h which
    // are identical to those in such processes as HydroMechanics,
    // RichardsMechanics, and etc, into a common place.
    ProcessLib::StokesFlow::createLocalAssemblers<GlobalDim,
                                                  LocalAssemblerData>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        pv.getShapeFunctionOrder(), _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

    _process_data.pressure_interpolated =
        MeshLib::getOrCreateMeshProperty<double>(
            const_cast<MeshLib::Mesh&>(mesh), "pressure_interpolated",
            MeshLib::MeshItemType::Node, 1);
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::initializeBoundaryConditions()
{
    assert(_use_monolithic_scheme);
    {
        int const process_id = 0;
        initializeProcessBoundaryConditionsAndSourceTerms(
            *_local_to_global_index_map, process_id);
    }
}

template <int GlobalDim>
MathLib::MatrixSpecifications
StokesFlowProcess<GlobalDim>::getMatrixSpecifications(
    const int /*process_id*/) const
{
    assert(_use_monolithic_scheme);
    {
        auto const& l = *_local_to_global_index_map;

        return {l.dofSizeWithoutGhosts(), l.dofSizeWithoutGhosts(),
                &l.getGhostIndices(), &this->_sparsity_pattern};
    }
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::assembleConcreteProcess(
    const double t, double const dt, std::vector<GlobalVector*> const& x,
    std::vector<GlobalVector*> const& xdot, int const process_id,
    GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b)
{
    DBUG("Assemble StokesFlowProcess.");

    std::vector<std::reference_wrapper<NumLib::LocalToGlobalIndexMap>>
        dof_tables;
    assert(_use_monolithic_scheme);
    {
        dof_tables.push_back(std::ref(*_local_to_global_index_map));
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    // Call global assembler for each local assembly item.
    GlobalExecutor::executeSelectedMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, t, dt, x, xdot, process_id, M, K,
        b);
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::assembleWithJacobianConcreteProcess(
    const double /*t*/, double const /*dt*/,
    std::vector<GlobalVector*> const& /*x*/,
    std::vector<GlobalVector*> const& /*xdot*/, const double /*dxdot_dx*/,
    const double /*dx_dx*/, int const /*process_id*/, GlobalMatrix& /*M*/,
    GlobalMatrix& /*K*/, GlobalVector& /*b*/, GlobalMatrix& /*Jac*/)
{
    OGS_FATAL(
        "Assembly of Jacobian matrix has not yet been implemented for "
        "StokesFlowProcess.");
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::computeSecondaryVariableConcrete(
    double const t,
    double const dt,
    std::vector<GlobalVector*> const& x,
    GlobalVector const& x_dot,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    assert(_use_monolithic_scheme);
    {
        dof_tables.push_back(_local_to_global_index_map.get());
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &StokesFlowLocalAssemblerInterface::computeSecondaryVariable,
        _local_assemblers, pv.getActiveElementIDs(), dof_tables, t, dt, x,
        x_dot, process_id);
}

template <int GlobalDim>
void StokesFlowProcess<GlobalDim>::postTimestepConcreteProcess(
    std::vector<GlobalVector*> const& x,
    const double t,
    const double dt,
    int const process_id)
{
    if (process_id != 0)
    {
        return;
    }

    std::vector<NumLib::LocalToGlobalIndexMap const*> dof_tables;
    dof_tables.reserve(x.size());
    assert(_use_monolithic_scheme);
    {
        dof_tables.push_back(_local_to_global_index_map.get());
    }

    ProcessLib::ProcessVariable const& pv = getProcessVariables(process_id)[0];
    GlobalExecutor::executeSelectedMemberOnDereferenced(
        &StokesFlowLocalAssemblerInterface::postTimestep, _local_assemblers,
        pv.getActiveElementIDs(), dof_tables, x, t, dt);
}

template <int GlobalDim>
NumLib::LocalToGlobalIndexMap const& StokesFlowProcess<GlobalDim>::getDOFTable(
    const int /*process_id*/) const
{
    assert(_use_monolithic_scheme);
    {
        return *_local_to_global_index_map;
    }
}

template class StokesFlowProcess<2>;
}  // namespace StokesFlow
}  // namespace ProcessLib
