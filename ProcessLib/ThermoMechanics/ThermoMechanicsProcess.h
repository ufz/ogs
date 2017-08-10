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

#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/Process.h"
#include "ProcessLib/SmallDeformation/CreateLocalAssemblers.h"

#include "ThermoMechanicsFEM.h"
#include "ThermoMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoMechanics
{
template <int DisplacementDim>
class ThermoMechanicsProcess final : public Process
{
    using Base = Process;

public:
    ThermoMechanicsProcess(
        MeshLib::Mesh& mesh,
        std::unique_ptr<ProcessLib::AbstractJacobianAssembler>&&
            jacobian_assembler,
        std::vector<std::unique_ptr<ParameterBase>> const& parameters,
        unsigned const integration_order,
        std::vector<std::reference_wrapper<ProcessVariable>>&&
            process_variables,
        ThermoMechanicsProcessData<DisplacementDim>&& process_data,
        SecondaryVariableCollection&& secondary_variables,
        NumLib::NamedFunctionCaller&& named_function_caller)
        : Process(mesh, std::move(jacobian_assembler), parameters,
                  integration_order, std::move(process_variables),
                  std::move(secondary_variables),
                  std::move(named_function_caller)),
          _process_data(std::move(process_data))
    {
    }

    //! \name ODESystem interface
    //! @{

    bool isLinear() const override { return false; }
    //! @}

private:
    using LocalAssemblerInterface = ThermoMechanicsLocalAssemblerInterface;

    void initializeConcreteProcess(
        NumLib::LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh,
        unsigned const integration_order) override
    {
        ProcessLib::SmallDeformation::createLocalAssemblers<
        DisplacementDim, ThermoMechanicsLocalAssembler>(
        mesh.getElements(), dof_table, _local_assemblers,
        mesh.isAxiallySymmetric(), integration_order, _process_data);

        // TODO move the two data members somewhere else.
        // for extrapolation of secondary variables
        std::vector<MeshLib::MeshSubsets> all_mesh_subsets_single_component;
        all_mesh_subsets_single_component.emplace_back(
            _mesh_subset_all_nodes.get());
        _local_to_global_index_map_single_component.reset(
            new NumLib::LocalToGlobalIndexMap(
                std::move(all_mesh_subsets_single_component),
                // by location order is needed for output
                NumLib::ComponentOrder::BY_LOCATION));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xx",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaXX));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaYY));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "sigma_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaXY));

        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "sigma_xz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaXZ));

            Base::_secondary_variables.addSecondaryVariable(
                "sigma_yz",
                makeExtrapolator(
                    1, getExtrapolator(), _local_assemblers,
                    &ThermoMechanicsLocalAssemblerInterface::getIntPtSigmaYZ));
        }
        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xx",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtEpsilonXX));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_yy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtEpsilonYY));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_zz",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtEpsilonZZ));

        Base::_secondary_variables.addSecondaryVariable(
            "epsilon_xy",
            makeExtrapolator(
                1, getExtrapolator(), _local_assemblers,
                &ThermoMechanicsLocalAssemblerInterface::getIntPtEpsilonXY));
        if (DisplacementDim == 3)
        {
            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_yz",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                                 &ThermoMechanicsLocalAssemblerInterface::
                                     getIntPtEpsilonYZ));

            Base::_secondary_variables.addSecondaryVariable(
                "epsilon_xz",
                makeExtrapolator(1, getExtrapolator(), _local_assemblers,
                                 &ThermoMechanicsLocalAssemblerInterface::
                                     getIntPtEpsilonXZ));
        }
    }

    void assembleConcreteProcess(
        const double t, GlobalVector const& x, GlobalMatrix& M, GlobalMatrix& K,
        GlobalVector& b, StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("Assemble ThermoMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assemble,
            _local_assemblers, *_local_to_global_index_map, t, x, M, K, b,
            coupling_term);
    }

    void assembleWithJacobianConcreteProcess(
        const double t, GlobalVector const& x, GlobalVector const& xdot,
        const double dxdot_dx, const double dx_dx, GlobalMatrix& M,
        GlobalMatrix& K, GlobalVector& b, GlobalMatrix& Jac,
        StaggeredCouplingTerm const& coupling_term) override
    {
        DBUG("AssembleJacobian ThermoMechanicsProcess.");

        // Call global assembler for each local assembly item.
        GlobalExecutor::executeMemberDereferenced(
            _global_assembler, &VectorMatrixAssembler::assembleWithJacobian,
            _local_assemblers, *_local_to_global_index_map, t, x, xdot,
            dxdot_dx, dx_dx, M, K, b, Jac, coupling_term);
    }

    void preTimestepConcreteProcess(GlobalVector const& x, double const t,
                                    double const dt) override
    {
        DBUG("PreTimestep ThermoMechanicsProcess.");

        _process_data.dt = dt;
        _process_data.t = t;

        GlobalExecutor::executeMemberOnDereferenced(
            &ThermoMechanicsLocalAssemblerInterface::preTimestep,
            _local_assemblers, *_local_to_global_index_map, x, t, dt);
    }

    void postTimestepConcreteProcess(GlobalVector const& x) override
    {
        DBUG("PostTimestep ThermoMechanicsProcess.");

        GlobalExecutor::executeMemberOnDereferenced(
            &ThermoMechanicsLocalAssemblerInterface::postTimestep,
            _local_assemblers, *_local_to_global_index_map, x);
    }

private:
    std::vector<MeshLib::Node*> _base_nodes;
    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_base_nodes;
    ThermoMechanicsProcessData<DisplacementDim> _process_data;

    std::vector<std::unique_ptr<LocalAssemblerInterface>> _local_assemblers;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap>
        _local_to_global_index_map_single_component;
};

}  // namespace ThermoMechanics
}  // namespace ProcessLib
