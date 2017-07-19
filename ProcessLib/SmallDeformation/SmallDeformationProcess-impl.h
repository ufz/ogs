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

#include "BaseLib/Functional.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "SmallDeformationFEM.h"

namespace ProcessLib
{
namespace SmallDeformation
{
template <int DisplacementDim>
SmallDeformationProcess<DisplacementDim>::SmallDeformationProcess(
    MeshLib::Mesh& mesh,
    std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&& jacobian_assembler,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters,
    unsigned const integration_order,
    std::vector<std::reference_wrapper<ProcessVariable>>&& process_variables,
    SmallDeformationProcessData<DisplacementDim>&& process_data,
    SecondaryVariableCollection&& secondary_variables,
    NumLib::NamedFunctionCaller&& named_function_caller)
    : Process(mesh, std::move(jacobian_assembler), parameters,
              integration_order, std::move(process_variables),
              std::move(secondary_variables), std::move(named_function_caller)),
      _process_data(std::move(process_data))
{
    _nodal_forces = MeshLib::getOrCreateMeshProperty<double>(
        mesh, "NodalForces", MeshLib::MeshItemType::Node, DisplacementDim);
}

template <int DisplacementDim>
bool SmallDeformationProcess<DisplacementDim>::isLinear() const
{
    return false;
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::initializeConcreteProcess(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh,
    unsigned const integration_order)
{
    ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, SmallDeformationLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

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
    _nodal_forces->resize(DisplacementDim * mesh.getNumberOfNodes());

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXX));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaYY));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "sigma_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtSigmaXY));

    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaXZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtSigmaYZ));
    }

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xx",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXX));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_yy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonYY));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_zz",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonZZ));

    Base::_secondary_variables.addSecondaryVariable(
        "epsilon_xy",
        makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                         &LocalAssemblerInterface::getIntPtEpsilonXY));
    if (DisplacementDim == 3)
    {
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonYZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xz",
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             &LocalAssemblerInterface::getIntPtEpsilonXZ));
    }

    auto const int_vars = _process_data.material->getInternalVariables();
    INFO("SmallDef has %lu internal vars.", int_vars.size());
    for (auto const& name_fct : int_vars)
    {
        auto const& name = name_fct.first;
        auto const& fct = name_fct.second;
        INFO("internal var %s.", name.c_str());

        auto getIntPtValues = BaseLib::easyBind(
            [fct](LocalAssemblerInterface const& loc_asm,
                  const double /*t*/,
                  GlobalVector const& /*current_solution*/,
                  NumLib::LocalToGlobalIndexMap const& /*dof_table*/,
                  std::vector<double>& cache) -> std::vector<double> const& {
                const unsigned num_int_pts =
                    loc_asm.getNumberOfIntegrationPoints();

                cache.clear();
                cache.reserve(num_int_pts);

                for (unsigned i = 0; i < num_int_pts; ++i)
                {
                    auto const& state = loc_asm.getMaterialStateVariablesAt(i);
                    cache.push_back(fct(state));
                }

                return cache;
            });

        Base::_secondary_variables.addSecondaryVariable(
            name,
            makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                             std::move(getIntPtValues)));
    }
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::assembleConcreteProcess(
    const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
    GlobalVector& b, StaggeredCouplingTerm const& coupling_term)
{
    DBUG("Assemble SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assemble, _local_assemblers,
        *_local_to_global_index_map, t, x, M, K, b, coupling_term);
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::
    assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term)
{
    DBUG("AssembleWithJacobian SmallDeformationProcess.");

    // Call global assembler for each local assembly item.
    GlobalExecutor::executeMemberDereferenced(
        _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
        _local_assemblers, *_local_to_global_index_map, t, x, xdot, dxdot_dx,
        dx_dx, M, K, b, Jac, coupling_term);

    b.copyValues(*_nodal_forces);
    std::transform(_nodal_forces->begin(), _nodal_forces->end(),
                   _nodal_forces->begin(), [](double val) { return -val; });
}

template <int DisplacementDim>
void SmallDeformationProcess<DisplacementDim>::preTimestepConcreteProcess(
    GlobalVector const& x, double const t, double const dt)
{
    DBUG("PreTimestep SmallDeformationProcess.");

    _process_data.dt = dt;
    _process_data.t = t;

    GlobalExecutor::executeMemberOnDereferenced(
        &LocalAssemblerInterface::preTimestep, _local_assemblers,
        *_local_to_global_index_map, x, t, dt);
}

}  // namespace SmallDeformation
}  // namespace ProcessLib
