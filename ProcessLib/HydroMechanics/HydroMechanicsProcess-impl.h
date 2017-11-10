/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cassert>

#include "MeshLib/Elements/Utils.h"
#include "ProcessLib/HydroMechanics/CreateLocalAssemblers.h"
#include "ProcessLib/Process.h"

#include "HydroMechanicsFEM.h"
#include "HydroMechanicsProcessData.h"

namespace ProcessLib
{
namespace HydroMechanics
{
template <int DisplacementDim>
HydroMechanicsProcess<DisplacementDim>::HydroMechanicsProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::vector<std::reference_wrapper<ProcessVariable>>>&&
        process_variables,
    HydroMechanicsProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
}

template <int DisplacementDim>
bool HydroMechanicsProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::constructDofTable()
{
    // Create single component dof in every of the mesh's nodes.
    _mesh_subset_all_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_mesh.getNodes());
    // Create single component dof in the mesh's base nodes.
    _base_nodes = MeshLib::getBaseNodes(_mesh.getElements());
    _mesh_subset_base_nodes =
        std::make_unique<MeshLib::MeshSubset>(_mesh, &_base_nodes);

    // Collect the mesh subsets in a vector.

    // For pressure, which is the first
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets;
    all_mesh_subsets.emplace_back(_mesh_subset_base_nodes.get());

    // For displacement.
    std::generate_n(
        std::back_inserter(all_mesh_subsets),
        getProcessVariables()[1].get().getNumberOfComponents(),
        [&]() { return MeshLib::MeshSubsets{_mesh_subset_all_nodes.get()}; });

    std::vector<int> const vec_n_components{1, DisplacementDim};
    _local_to_global_index_map =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets), vec_n_components,
            NumLib::ComponentOrder::BY_LOCATION);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::HydroMechanics::createLocalAssemblers<
        DisplacementDim, HydroMechanicsLocalAssembler>(
        mesh.getDimension(), mesh.getElements(), dof_table,
        // use displacment process variable for shapefunction order
        getProcessVariables()[1].get().getShapeFunctionOrder(),
        _local_assemblers, mesh.isAxiallySymmetric(), integration_order,
        _process_data);

    // TODO move the two data members somewhere else.
    // for extrapolation of secondary variables
    std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
    all_mesh_subsets_single_component.emplace_back(
        _mesh_subset_all_nodes.get());
    _local_to_global_index_map_single_component =
        std::make_unique<NumLib::LocalToGlobalIndexMap>(
            std::move(all_mesh_subsets_single_component),
            // by location order is needed for output
            NumLib::ComponentOrder::BY_LOCATION);

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &LocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(
            1, getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtEpsilonXY));

    Base::_secondary_variables.addSecondaryVariable(
        "velocity",
        makeExtrapolator(
            mesh.getDimension(), getExtrapolator(), _local_assemblers,
            &LocalAssemblerInterface::getIntPtDarcyVelocity));
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b)
{
    DBUG("Assemble HydroMechanicsProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, _coupled_solutions);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac)
{
    DBUG("AssembleJacobian HydroMechanicsProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, _coupled_solutions);
}
template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt,
    const int /*process_id*/)
{
    DBUG("PreTimestep HydroMechanicsProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

template <int DisplacementDim>
void HydroMechanicsProcess<DisplacementDim>::postTimestepConcreteProcess(
    GlobalVector const& x)
{
    DBUG("PostTimestep HydroMechanicsProcess.");

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::postTimestep, _local_assemblers,
        *_local_to_global_index_map, x);
}

}  // namespace HydroMechanics
}  // namespace ProcessLib
